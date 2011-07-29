#include <cassert>
#include <algorithm>

#include "sam_buffer.h"
#include "sam_helper.h"
#include "sam_order.h"
#include "readsim_aux.h"

/*
  The sam_order instance of 'less' or 'equal' caches which function to
  use.  The calling instance (here the same instance) caches the
  contig_offsets / contig_lengths for use by the dereferenced function
  pointer.
 */

//sort on whichever 'less' function is chosen by 'sam_order'
bool LessSAMLinePtr::operator()(SamLine const* a,
                                SamLine const* b)
{
    return (this->sam_order->*(this->sam_order->less))(*a, *b);
}

//sorts paired sam records with the same logic as the above
//'LessSAMLinePtr' hierarchically
bool LessSAMLinePair::operator()(SAMPTR_PAIR const& a,
                                 SAMPTR_PAIR const& b)
{
    return (this->sam_order->*(this->sam_order->less))(*a.first, *b.first)
        || ((this->sam_order->*(this->sam_order->equal))(*a.first, *b.first)
            && (this->sam_order->*(this->sam_order->less))(*a.second, *b.second));
}


//sort on fragment ID (qid or qname)
bool LessSAMLineFragmentIDPtr::operator()(SamLine const* a,
                                          SamLine const* b)
{
    return this->sam_order->less_fid(*a, *b);
}


//sort on flattened position
bool LessSAMLinePtrMatePair::operator()(SamLine const* a,
                                        SamLine const* b)
{
    return this->sam_order->less_position_rid(*a, *b);
}


SamBuffer::SamBuffer(SamOrder const* _sam_order, bool pairs_same_strand) : 
    
    output_pairs_as_same_strand(pairs_same_strand),
    less_ordering(_sam_order),
    less_entry_matepair(_sam_order),
    sam_order(_sam_order),
    low_bound(NULL),
    unique_entry_pairs(LessSAMLinePair(_sam_order)),
    yet_unpaired_entries(LessSAMLineFragmentIDPtr(_sam_order)),
    unpaired_entries(LessSAMLinePtr(_sam_order))
{
}



//advance lowbound with proposed new lowbound, if it is indeed
//greater.  if not, do nothing and return false.
bool SamBuffer::safe_advance_lowbound(SamLine const* _proposed_new_low_bound)
{
    if (this->low_bound == NULL
        || this->less_ordering(this->low_bound, _proposed_new_low_bound))
    {
        this->low_bound = _proposed_new_low_bound;
        return true;
    }
    else
    {
        return false;
    }
}

void SamBuffer::update_lowbound(SamLine const* _low_bound) 
{ 
    //the low bound should never be decreased
    assert(this->low_bound == NULL
           || this->less_ordering(this->low_bound, _low_bound));

    this->low_bound = _low_bound; 
}



//returns pair, with a pointer to the inserted entry or a pre-existing
//duplicate entry, and a bool indicating whether insertion was
//successful. This is the same logic as with STL container insert()
//but with a pointer to the actual item rather than an iterator
std::pair<SamLine const*, bool>
SamBuffer::insert(SamLine const* entry)
{
    std::pair<SamLine const*, bool> result(static_cast<SamLine const*>(NULL), false);

    if (! entry->paired_in_sequencing())
    {
        std::pair<SINGLE_READ_SET::iterator, bool> insert = 
            this->unpaired_entries.insert(entry);
        assert(insert.second);
        result = std::make_pair(entry, true);
    }
    else
    {
        //is mate stored already?  

        //will find all as-yet unpaired entries with the same qname as
        //the entry in the regime where the first and then second of
        //the entries are inserted, this will trigger a non-empty range
        //every time the second in the pair is inserted.
        std::pair<SINGLE_READ_ID_MSET::const_iterator,
                  SINGLE_READ_ID_MSET::const_iterator> yet_unpaired_range =
            this->yet_unpaired_entries.equal_range(entry);

        //it is an error if the existing entry finds more than one mapped mate
        bool mate_found = false;
        //bool duplicate_pair_found = false;

        //
        // used to determine whether the current entry or existing
        // unpaired entry should be first in the pair.

        SINGLE_READ_ID_MSET::iterator unpaired_iter = yet_unpaired_range.first;
        SINGLE_READ_ID_MSET::iterator unpaired_next_iter = unpaired_iter;

        //go through each unpaired entry that is a possible match.
        //policy:  once a match is found, record it and break from the loop.
        //any other possible matches will be searched by subsequent suitors
        while (unpaired_iter != yet_unpaired_range.second)
        {
            SamLine const* yet_unpaired_entry = *unpaired_iter;
            unpaired_next_iter = unpaired_iter;
            ++unpaired_next_iter;

            bool is_mapped_mate_pair = AreMappedMatePairs(*entry, *yet_unpaired_entry);
            bool is_unmapped_mate_pair = AreUnmappedMatePairs(*entry, *yet_unpaired_entry);
            if (is_mapped_mate_pair || is_unmapped_mate_pair)
            {
                mate_found = true;

                bool unpaired_iter_is_leftmost = this->less_entry_matepair(yet_unpaired_entry, entry);
                SamLine const* left = unpaired_iter_is_leftmost ? yet_unpaired_entry : entry;
                SamLine const* right = unpaired_iter_is_leftmost ? entry : yet_unpaired_entry;
                assert(left != right);
                
                //now, in accordance with sam format spec
                if (left->isize < 0 || right->isize > 0)
                {
                    //bool unpaired_iter_is_leftmost = this->less_entry_matepair(yet_unpaired_entry, entry);
                    fprintf(stderr, "SAM format error: inappropriate 'isize' fields.  Should be "
                            "positive for left and negative for right-most reads\n"
                            "Entries:\n");
                    left->print(stderr, true);
                    right->print(stderr, true);
                    exit(1);
                }
                
                std::pair<PAIRED_READ_SET::const_iterator, bool> insert_success = 
                    this->unique_entry_pairs.insert(std::make_pair(left, right));
                
                if (insert_success.second)
                {
                    //successfully inserted the pair together.
                    result = std::make_pair(entry, true);
                    this->yet_unpaired_entries.erase(unpaired_iter);
                }
                else
                {
                    //this entry reveals that it and its mate create a
                    //duplicate read pair in the set of
                    //unique_entry_pairs.
                    SamLine const* pre_existing_entry = 
                        unpaired_iter_is_leftmost 
                        ? (*insert_success.first).second
                        : (*insert_success.first).first;

                    result = std::make_pair(pre_existing_entry, false);

                    if (this->low_bound == yet_unpaired_entry)
                    {
                        //we want to then set the low bound to the SamLine
                        //that represents the pre-existing duplicate SamLine
                        //that the low_bound was previously set to.
                        this->low_bound = unpaired_iter_is_leftmost ?
                            (*insert_success.first).first :
                            (*insert_success.first).second;
                        assert(! this->less_entry_matepair(this->low_bound, yet_unpaired_entry));
                        assert(! this->less_entry_matepair(yet_unpaired_entry, this->low_bound));
                        
                    }
                    this->yet_unpaired_entries.erase(unpaired_iter);
                    delete yet_unpaired_entry;
                    assert(entry != this->low_bound);
                    delete entry;
                }
                break;
            }
            unpaired_iter = unpaired_next_iter;
        }
        if (! mate_found)
        {
            //treat as yet unpaired
            this->yet_unpaired_entries.insert(entry);
            result = std::make_pair(entry, true);
        }
    }
    assert(result.first != NULL);

    return result;
}

/*
  print entries preceding low_bound, and delete them from memory.  if
  any of the filehandles are NULL, do not print to that file, but do
  everything else.

  sets low_bound to NULL after the purge.
 */
void SamBuffer::purge(FILE * output_sam_fh, 
                      FILE * output_first_fastq_fh, 
                      FILE * output_second_fastq_fh, 
                      bool ignore_bound)
{

    //since no new read can come before low-bound, we can be sure that
    //all completed pairs in which both reads are before the low bound
    //can be purged.

    //what about the yet_unpaired entries?  how do we know when to
    //purge them?  we can only know to purge a yet-unpaired entry if
    //we know that its potential mate has to be below the low bound
    
    //to do this, we need an ordering that keeps each pair of entries
    //together, and sorts them by the overall fragment ordering.

    //what would this be?  it would have to take the minimum alignment
    //position between it and the mate.

    PAIRED_READ_SET::const_iterator paired_bound;
    SINGLE_READ_SET::const_iterator single_bound;
    SINGLE_READ_SET::const_iterator unpaired_bound;
    //SINGLE_READ_ID_MSET::const_iterator yet_unpaired_bound;

    SamLine const* temp_low_bound = this->low_bound;

    if (ignore_bound)
    {
        paired_bound = this->unique_entry_pairs.end();
        // single_bound = this->single_entries_from_pairs.end();
        unpaired_bound = this->unpaired_entries.end();
        //yet_unpaired_bound = this->yet_unpaired_entries.end();

        this->low_bound = NULL;
    }

    else
    {
        //pay attention to the bound that is set.  why wouldn't we
        //want to purge in the same order as the SamBuffer was
        //specified?
        if (this->low_bound == NULL)
        {
            //do not purge anything if low_bound is NULL
            paired_bound = this->unique_entry_pairs.begin();
            unpaired_bound = this->unpaired_entries.begin();
            //yet_unpaired_bound = this->yet_unpaired_entries.begin();
        }
        else
        {
            SamLine const* min_incomplete =
                this->yet_unpaired_entries.empty()
                ? this->low_bound
                : std::min(this->low_bound, *this->yet_unpaired_entries.begin(), this->less_ordering);

            paired_bound =
                this->unique_entry_pairs.lower_bound
                (std::make_pair(min_incomplete, min_incomplete));

            // if (paired_bound != this->unique_entry_pairs.begin())
            // {
            //     //this allows us to recover an existing SAMline pair
            //     //entry and thus avoid the following problem: Given a
            //     //pre-existing set of pairs (A B) (C D), the
            //     //pseudo-lower-bound pair (B B) would insert as: (A B)
            //     //(B B) (C D). So, if we purge, we would delete (A B)
            //     //and thus delete the individual entry B, invalidating
            //     //the low_bound of (B B)
            //     --paired_bound;
            // }

            unpaired_bound =
                this->unpaired_entries.lower_bound(this->low_bound);

            //yet_unpaired_bound = 
            //this->yet_unpaired_entries.lower_bound(this->low_bound);

        }
    }

    if (output_first_fastq_fh != NULL)
    {
        FILE * output_second_fastq_used_fh =
            output_second_fastq_fh == NULL ? output_first_fastq_fh :
            output_second_fastq_fh;

        PAIRED_READ_SET::const_iterator pair_iter;
        for (pair_iter = this->unique_entry_pairs.begin();
             pair_iter != paired_bound; ++pair_iter)
        {
            print_paired_fastq_entries(output_first_fastq_fh,
                                       output_second_fastq_used_fh,
                                       (*pair_iter).first,
                                       (*pair_iter).second);
        }
    }


    PAIRED_READ_SET::const_iterator pair_iter;

    for (pair_iter = this->unique_entry_pairs.begin();
         pair_iter != paired_bound; ++pair_iter)
    {
        SamLine const* first = (*pair_iter).first;
        SamLine const* second = (*pair_iter).second;

        if (output_sam_fh != NULL)
        {
            first->print(output_sam_fh, 
                         this->output_pairs_as_same_strand 
                         && first->second_read_in_pair());

            second->print(output_sam_fh, 
                          this->output_pairs_as_same_strand 
                          && second->second_read_in_pair());
        }

        //since the paired_bound contains the
        // assert(first != this->low_bound);
        // assert(second != this->low_bound);

        delete first;
        delete second;
    }

    //fprintf(stderr, "unique_read_pairs.size(): %Zu ...", this->unique_read_pairs.size());
    this->unique_entry_pairs.erase(this->unique_entry_pairs.begin(), paired_bound);
    //fprintf(stderr, "%Zu\n", this->unique_entry_pairs.size());
    //fflush(stderr);

    SINGLE_READ_SET::const_iterator read_iter;

    for (read_iter = this->unpaired_entries.begin();
         read_iter != unpaired_bound; ++read_iter)
    {
        SamLine const* cur_read = *read_iter;
        bool flip_query_strand_flag = this->output_pairs_as_same_strand &&
            cur_read->second_read_in_pair();

        if (output_sam_fh != NULL)
        {
            cur_read->print(output_sam_fh, flip_query_strand_flag);
        }
        
        // assert(*read_iter != this->low_bound);
        assert(this->yet_unpaired_entries.find(cur_read) == this->yet_unpaired_entries.end());
        delete cur_read;
    }
    //fprintf(stderr, "unpaired_entries.size(): %Zu ...", this->unpaired_entries.size());
    this->unpaired_entries.erase(this->unpaired_entries.begin(), unpaired_bound);
    //fprintf(stderr, "%Zu\n", this->unpaired_entries.size());
    //fflush(stderr);


    //output all yet_unpaired_entries that fall below the bound.  why?
    /*
    SINGLE_READ_ID_MSET::iterator unpaired_iter, iter_next;
    iter_next = this->yet_unpaired_entries.begin();

    while (iter_next != yet_unpaired_bound)
    {

        unpaired_iter = iter_next;
        SamLine * cur_read = const_cast<SamLine *>(*unpaired_iter);
        ++iter_next;

        bool flip_query_strand_flag = this->output_pairs_as_same_strand &&
            cur_read->second_read_in_pair();
        
        if (ignore_bound || this->less_entry_matepair(cur_read, this->low_bound))
        {
            cur_read->flag ^= SamFlags::MAPPED_IN_PROPER_PAIR;
            if (output_sam_fh != NULL)
            {
                cur_read->print(output_sam_fh, flip_query_strand_flag);
            }
            this->yet_unpaired_entries.erase(unpaired_iter);
            delete cur_read;
        }
    }
    */

    if (ignore_bound)
    {
        //restore old low_bound
        this->low_bound = temp_low_bound;

        if (! this->yet_unpaired_entries.empty())
        {
            fprintf(stderr, "Warning: SamBuffer::purge(): request to purge all remaining "
                    "entries yet there are still %Zu as-yet unpaired entries\n",
                    this->yet_unpaired_entries.size());
            exit(1);
        }
    }
    
    this->low_bound = NULL;

    if (output_sam_fh != NULL)
    {
        fflush(output_sam_fh);
    }
    if (output_first_fastq_fh != NULL)
    {
        fflush(output_first_fastq_fh);
    }
    if (output_second_fastq_fh != NULL)
    {
        fflush(output_second_fastq_fh);
    }

}
