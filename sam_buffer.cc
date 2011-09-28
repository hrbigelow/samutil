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
                                SamLine const* b) const
{
    return (this->sam_order->*(this->sam_order->less))(*a, *b);
}


//sort on fragment ID (qid or qname)
bool LessSAMLineFragmentIDPtr::operator()(SamLine const* a,
                                          SamLine const* b) const
{
    return this->sam_order->less_fragment_id(*a, *b);
}


//sort on flattened position
bool LessSAMLinePtrMatePair::operator()(SamLine const* a,
                                        SamLine const* b) const
{
    return this->sam_order->less_position_rid(*a, *b);
}


SamBuffer::SamBuffer(SamOrder const* _sam_order, bool pairs_same_strand) : 
    
    output_pairs_as_same_strand(pairs_same_strand),
    less_ordering(_sam_order),
    less_entry_matepair(_sam_order),
    sam_order(_sam_order),
    unique_entries(LessSAMLinePtr(_sam_order)),
    incomplete_entries(LessSAMLineFragmentIDPtr(_sam_order))
{
}



//advance lowbound with proposed new lowbound, if it is indeed
//greater.  if not, do nothing and return false.
/*
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
*/

InsertResult::InsertResult(SINGLE_READ_SET::iterator _se, bool _wi, SamLine * _sam)
    : surviving_entry(_se), was_inserted(_wi), remaining_entry(_sam) { }

InsertResult::InsertResult() : 
    surviving_entry(), was_inserted(false), remaining_entry(NULL) { }


//returns pair, with a pointer to the inserted entry or a pre-existing
//duplicate entry, and a bool indicating whether insertion was
//successful. This is the same logic as with STL container insert()
//but with a pointer to the actual item rather than an iterator
InsertResult SamBuffer::insert(SamLine * entry)
{
    InsertResult result;

    if (entry->flag.is_rsam_format)
    {
        //entry is already in rSAM format.  Simply insert
        std::pair<SINGLE_READ_SET::iterator, bool> ins = 
            this->unique_entries.insert(entry);

        result = InsertResult(ins.first, ins.second, entry);
    }
    else if ((! entry->flag.is_rsam_format) && ! entry->flag.multi_fragment_template)
    {
        //we have traditional SAM, and it is a singleton.  convert it to rSAM
        SamLine const* line_array[] = { entry };
        SamLine * merged = new SamLine(line_array, 1, SamLine::expected_read_layout);
        delete entry;

        //this entry is 
        std::pair<SINGLE_READ_SET::iterator, bool> ins = 
            this->unique_entries.insert(merged);

        result = InsertResult(ins.first, ins.second, merged);
    }
    else
    {
        //traditional SAM format, multi-fragment template entry.  Needs to be
        //joined with its mate and merged.

        //is mate stored already?  

        //will find all as-yet unpaired entries with the same qname as
        //the entry in the regime where the first and then second of
        //the entries are inserted, this will trigger a non-empty range
        //every time the second in the pair is inserted.
        std::pair<SINGLE_READ_ID_MSET::iterator,
                  SINGLE_READ_ID_MSET::iterator> yet_unmerged_range =
            this->incomplete_entries.equal_range(entry);

        //it is an error if the existing entry finds more than one mapped mate
        bool mate_found = false;
        //bool duplicate_pair_found = false;

        // used to determine whether the current entry or existing
        // unmerged entry should be first in the pair.

        SINGLE_READ_ID_MSET::iterator unmerged_iter = yet_unmerged_range.first;
        SINGLE_READ_ID_MSET::iterator unmerged_next_iter = unmerged_iter;

        //go through each unmerged entry that is a possible match.
        //policy:  once a match is found, record it and break from the loop.
        //any other possible matches will be searched by subsequent suitors
        while (unmerged_iter != yet_unmerged_range.second)
        {
            SamLine const* yet_unmerged_entry = *unmerged_iter;
            unmerged_next_iter = unmerged_iter;
            ++unmerged_next_iter;

            bool is_mapped_mate_pair = AreMappedMatePairs(*entry, *yet_unmerged_entry);
            bool is_unmapped_mate_pair = AreUnmappedMatePairs(*entry, *yet_unmerged_entry);
            if (is_mapped_mate_pair || is_unmapped_mate_pair)
            {
                mate_found = true;

                bool unmerged_iter_is_leftmost = this->less_entry_matepair(yet_unmerged_entry, entry);
                SamLine const* left = unmerged_iter_is_leftmost ? yet_unmerged_entry : entry;
                SamLine const* right = unmerged_iter_is_leftmost ? entry : yet_unmerged_entry;
                assert(left != right);
                
                //now, in accordance with sam format spec
                if (left->tlen < 0 || right->tlen > 0)
                {
                    //bool unmerged_iter_is_leftmost = this->less_entry_matepair(yet_unmerged_entry, entry);
                    fprintf(stderr, "SAM format error: inappropriate 'tlen' fields.  Should be "
                            "positive for left and negative for right-most reads\n"
                            "Entries:\n");
                    left->fprint(stderr);
                    right->fprint(stderr);
                    exit(1);
                }

                SamLine const* line_array[] = { left, right };

                SamLine * merged = new SamLine(line_array, 2, SamLine::expected_read_layout);

                delete left;
                delete right;
                this->incomplete_entries.erase(unmerged_iter);
                    
                std::pair<SINGLE_READ_SET::iterator, bool> 
                    ins = this->unique_entries.insert(merged);

                result = InsertResult(ins.first, ins.second, merged);

                break;
            }
            unmerged_iter = unmerged_next_iter;
        }
        if (! mate_found)
        {
            // treat as yet unmerged.  insertion will always succeed
            // because it is a multi-set.
            this->incomplete_entries.insert(entry);

            // how to we signal that the iterator is meaningless?
            // we don't.  just note that, if insertion succeeded, we don't expect
            // the caller to fuss with the iterator.
            result = InsertResult(this->unique_entries.end(), true, entry);
        }
    }

    return result;
}



void SamBuffer::replace(SINGLE_READ_SET::iterator iter,
                        SamLine * replacement_entry)
{
    SINGLE_READ_SET::iterator prec = iter;

    // correct for a valid prec iterator
    if (iter == this->unique_entries.begin())
    {
        prec = this->unique_entries.end();
    }
    else
    {
        --prec;
    }
    
    delete *iter;
    this->unique_entries.erase(iter);
    this->unique_entries.insert(prec, replacement_entry);
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
                      bool output_rsam_format,
                      SamLine const* low_bound)
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

    SINGLE_READ_SET::iterator unique_bound;
    SINGLE_READ_ID_MSET::iterator incomplete_bound;

    // SamLine const* temp_low_bound = this->low_bound;

    if (low_bound == NULL)
    {
        unique_bound = this->unique_entries.end();
        incomplete_bound = this->incomplete_entries.end();
    }

    else
    {
        //pay attention to the bound that is set.  why wouldn't we
        //want to purge in the same order as the SamBuffer was
        //specified?
        SamLine const* min_incomplete =
            this->incomplete_entries.empty()
            ? low_bound
            : std::min(low_bound, *this->incomplete_entries.begin(), this->less_ordering);
        
        unique_bound = this->unique_entries.lower_bound(min_incomplete);
        
        // if (unique_bound != this->unique_entries.begin())
        // {
        //     //this allows us to recover an existing SAMline pair
        //     //entry and thus avoid the following problem: Given a
        //     //pre-existing set of pairs (A B) (C D), the
        //     //pseudo-lower-bound pair (B B) would insert as: (A B)
        //     //(B B) (C D). So, if we purge, we would delete (A B)
        //     //and thus delete the individual entry B, invalidating
        //     //the low_bound of (B B)
        //     --unique_bound;
        // }

        incomplete_bound = this->incomplete_entries.lower_bound(low_bound);
        
    }

    if (output_first_fastq_fh != NULL)
    {
        FILE * output_second_fastq_used_fh =
            output_second_fastq_fh == NULL ? output_first_fastq_fh :
            output_second_fastq_fh;

        SINGLE_READ_SET::const_iterator iter;
        for (iter = this->unique_entries.begin(); iter != unique_bound; ++iter)
        {
            assert(false);
            // print_paired_fastq_entries(output_first_fastq_fh,
            //                            output_second_fastq_used_fh,
            //                            (*iter).first,
            //                            (*iter).second);
        }
    }


    SINGLE_READ_SET::const_iterator iter;

    for (iter = this->unique_entries.begin(); iter != unique_bound; ++iter)
    {
        SamLine const* samline = *iter;

        if (output_sam_fh != NULL)
        {
            if (output_rsam_format)
            {
                samline->fprint(output_sam_fh);
            }
            else
            {
                assert(false);
                //samline->print_sam(output_sam_fh);
            }
            // first->print(output_sam_fh, 
            //              this->output_pairs_as_same_strand 
            //              && first->second_read_in_pair());
        }
        delete samline;
    }

    //fprintf(stderr, "unique_read_pairs.size(): %Zu ...", this->unique_read_pairs.size());
    this->unique_entries.erase(this->unique_entries.begin(), unique_bound);
    //fprintf(stderr, "%Zu\n", this->unique_entries.size());
    //fflush(stderr);

    SINGLE_READ_ID_MSET::const_iterator read_iter;

    for (read_iter = this->incomplete_entries.begin();
         read_iter != incomplete_bound; ++read_iter)
    {
        SamLine const* samline = *read_iter;
        // bool flip_query_strand_flag = this->output_pairs_as_same_strand &&
        //     samline->second_read_in_pair();

        if (output_sam_fh != NULL)
        {
            if (output_rsam_format)
            {
                samline->fprint(output_sam_fh);
            }
            else
            {
                assert(false);
                //samline->print_sam(output_sam_fh);
            }
        }
        delete samline;
    }
    //fprintf(stderr, "incomplete_entries.size(): %Zu ...", this->incomplete_entries.size());
    this->incomplete_entries.erase(this->incomplete_entries.begin(), incomplete_bound);
    //fprintf(stderr, "%Zu\n", this->incomplete_entries.size());
    //fflush(stderr);


    if (low_bound == NULL)
    {
        if (! this->incomplete_entries.empty())
        {
            fprintf(stderr, "Warning: SamBuffer::purge(): request to purge all remaining "
                    "entries yet there are still %Zu as-yet unpaired entries\n",
                    this->incomplete_entries.size());
            exit(1);
        }
    }
    
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
