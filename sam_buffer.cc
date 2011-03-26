#include <cassert>
#include <algorithm>

#include "sam_buffer.h"
#include "sam_helper.h"
#include "readsim_aux.h"

//sort on [rname, pos, strand, cigar, and optionally qname]
bool LessSAMLinePtr::operator()(SamLine const* a,
                                SamLine const* b)
{
    switch (this->sam_order)
    {
    case SAM_RID: return a->less_rid(*b); break;
    case SAM_POSITION_RID: return a->less_position_rid(*b); break;
    case SAM_RID_POSITION: return a->less_rid_position(*b); break;
    case SAM_POSITION: return a->less_position(*b); break;
    case SAM_FID_POSITION: return a->less_fid_position(*b); break;
    default: 
        fprintf(stderr, "LessSAMLinePtr: unknown sort order\n");
        exit(1);
        break;
    }
}


//sort on [qname, pair]
bool LessSAMLineFragmentIDPtr::operator()(SamLine const* a,
                                          SamLine const* b)
{
    return a->less_fid(*b);
}


bool LessSAMLinePair::operator()(SAMPTR_PAIR const& a,
                                 SAMPTR_PAIR const& b)
{
    switch (this->sam_order)
    {
    case SAM_RID: 
        return 
            a.first->less_rid(*b.first)
            || (a.first->equal_rid(*b.first)
                && a.second->less_rid(*b.second));
        break;
    case SAM_POSITION_RID: 
        return 
            a.first->less_position_rid(*b.first)
            || (a.first->equal_position_rid(*b.first)
                && a.second->less_position_rid(*b.second));
        break;
    case SAM_RID_POSITION: 
        return 
            a.first->less_rid_position(*b.first)
            || (a.first->equal_rid_position(*b.first)
                && a.second->less_rid_position(*b.second));
        break;
    case SAM_POSITION:
        return
            a.first->less_position(*b.first)
            || (a.first->equal_position(*b.first)
                && a.second->less_position(*b.second));
        break;
    case SAM_FID_POSITION:
        return
            a.first->less_fid_position(*b.first)
            || (a.first->equal_fid_position(*b.first)
                && a.second->less_fid_position(*b.second));
        break;
    default: 
        fprintf(stderr, "LessSAMLinePtr: unknown sort order\n");
        exit(1);
        break;

    }
}


SamBuffer::SamBuffer(SAM_ORDER _sam_order,
                     bool pairs_same_strand,
                     bool output_ones_based,
                     bool ignore_dup_maps) : 
    
    sam_order(_sam_order),
    output_pairs_as_same_strand(pairs_same_strand),
    output_is_ones_based(output_ones_based),
    ignore_duplicate_mapped_pairs(ignore_dup_maps),
    low_bound(NULL),
    less_ordering(sam_order),
    unique_entry_pairs(LessSAMLinePair(sam_order)),
    unpaired_entries(LessSAMLinePtr(sam_order))
{ }



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



//returns whether successfully stored.  if not stored, deletes entry
void SamBuffer::insert(SamLine const* entry)
{
    if (! entry->paired_in_sequencing())
    {
        std::pair<SINGLE_READ_SET::iterator, bool> insert = 
            this->unpaired_entries.insert(entry);
        
        assert(insert.second);
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
        bool duplicate_pair_found = false;

        //
        // used to determine whether the current entry or existing
        // unpaired entry should be first in the pair.
        LessSAMLinePtr less_entry_matepair(SAM_FID_POSITION);

        SINGLE_READ_ID_MSET::iterator unpaired_iter = yet_unpaired_range.first;
        SINGLE_READ_ID_MSET::iterator unpaired_next_iter = unpaired_iter;

        //go through each unpaired entry that is a possible match
        while (unpaired_iter != yet_unpaired_range.second)
        {
            SamLine const* yet_unpaired_entry = *unpaired_iter;
            unpaired_next_iter = unpaired_iter;
            ++unpaired_next_iter;

            bool is_mapped_mate_pair = AreMappedMatePairs(*entry, *yet_unpaired_entry);
            bool is_unmapped_mate_pair = AreUnmappedMatePairs(*entry, *yet_unpaired_entry);
            if (is_mapped_mate_pair || is_unmapped_mate_pair)
            {
                if (mate_found)
                {
                    if (is_mapped_mate_pair)
                    {
                        if (this->ignore_duplicate_mapped_pairs)
                        {
                            this->yet_unpaired_entries.erase(unpaired_iter);
                            duplicate_pair_found = true;
                            //assume the duplicate mapped pair is
                            //properly paired, and discard both
                            //without trying to find another mate
                            break;
                        }
                        else
                        {
                            fprintf(stderr, 
                                    "Error: SamBuffer::insert: found duplicate "
                                    "mapped mate pairs for read %s\n", entry->qname);
                            exit(1);
                        }
                    }
                    else
                    {
                        //found a duplicate unmapped mate pair.
                        //assume that it is simply a duplicate because
                        //there is insufficient information to
                        //distinguish its exact mate from the mate of
                        //another.  pretend you didn't see it, so that
                        //it is left there for another unmapped pair
                        //to pair up with
                    }
                }
                else
                {
                    //mate not yet found.  continue searching
                    bool unpaired_iter_is_leftmost = less_entry_matepair(yet_unpaired_entry, entry);
                    SamLine const* left = unpaired_iter_is_leftmost ? yet_unpaired_entry : entry;
                    SamLine const* right = unpaired_iter_is_leftmost ? entry : yet_unpaired_entry;
                    assert(left != right);
                
                    //now, in accordance with sam format spec
                    if (left->isize < 0 || right->isize > 0)
                    {
                        bool unpaired_iter_is_leftmost = less_entry_matepair(yet_unpaired_entry, entry);
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
                        //successfully inserted the pair together.  do not test further
                        //because we want to allow 
                        mate_found = true;
                        this->yet_unpaired_entries.erase(unpaired_iter);
                        break;
                    }
                    else
                    {
                        //this entry reveals that it and its mate create a
                        //duplicate read pair.
                        if (this->low_bound == yet_unpaired_entry)
                        {
                            //we want to then set the low bound to the SamLine
                            //that represents the pre-existing duplicate SamLine
                            //that the low_bound was previously set to.
                            this->low_bound = unpaired_iter_is_leftmost ?
                                (*insert_success.first).first :
                                (*insert_success.first).second;
                            assert(! less_entry_matepair(this->low_bound, yet_unpaired_entry));
                            assert(! less_entry_matepair(yet_unpaired_entry, this->low_bound));

                        }
                        delete yet_unpaired_entry;
                        this->yet_unpaired_entries.erase(unpaired_iter);
                        duplicate_pair_found = true;
                        mate_found = false;
                    }                    
                }
            }
            unpaired_iter = unpaired_next_iter;
        }


        if (duplicate_pair_found)
        {
            assert(entry != this->low_bound);
            delete entry;
        }
        else
        {
            if (! mate_found)
            {
                //treat as yet unpaired
                this->yet_unpaired_entries.insert(entry);
            }
        }
    }
}

/*
  print entries preceding low_bound, and delete them from memory.
  if any of the filehandles are NULL, do not print to that file, but
  do everything else.
 */
void SamBuffer::purge(FILE * output_sam_fh, 
                      FILE * output_first_fastq_fh, 
                      FILE * output_second_fastq_fh, 
                      bool ignore_bound)
{

    //since no new read can come before low-bound, no new pair can
    //contain any read coming before low_bound.  therefore we can
    //delete all pairs that contain such entries.
    PAIRED_READ_SET::const_iterator paired_bound;
    SINGLE_READ_SET::const_iterator single_bound;
    SINGLE_READ_SET::const_iterator unpaired_bound;

    SamLine const* temp_low_bound = this->low_bound;

    if (ignore_bound)
    {
        //at the end of input, consider the bound to be infinite
        paired_bound = this->unique_entry_pairs.end();
        // single_bound = this->single_entries_from_pairs.end();
        unpaired_bound = this->unpaired_entries.end();
        this->low_bound = NULL;
    }

    else
    {
        //pay attention to the bound that is set.  why wouldn't we want to purge in
        //the same order as the SamBuffer was specified?

        //bool use_ids = true;

        paired_bound =
            this->unique_entry_pairs.lower_bound
            (std::make_pair(this->low_bound, this->low_bound));

        unpaired_bound =
            std::lower_bound(this->unpaired_entries.begin(), 
                             this->unpaired_entries.end(),
                             this->low_bound);
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
        // pair_iter_copy = pair_iter;
        // ++pair_iter; // so this will be valid!

        if (output_sam_fh != NULL)
        {
            first->print(output_sam_fh, 
                         this->output_pairs_as_same_strand 
                         && first->second_read_in_pair());

            second->print(output_sam_fh, 
                          this->output_pairs_as_same_strand 
                          && second->second_read_in_pair());
        }
        
        delete first;
        delete second;
        // this->unique_entry_pairs.erase(pair_iter_copy);
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
        
        assert(*read_iter != this->low_bound);
        delete cur_read;
    }
    //fprintf(stderr, "unpaired_entries.size(): %Zu ...", this->unpaired_entries.size());
    this->unpaired_entries.erase(this->unpaired_entries.begin(), unpaired_bound);
    //fprintf(stderr, "%Zu\n", this->unpaired_entries.size());
    //fflush(stderr);


    //output all yet_unpaired_entries that fall below the bound
    SINGLE_READ_ID_MSET::iterator unpaired_iter, iter_next;
    LessSAMLinePtr less_sam(this->sam_order);

    iter_next = this->yet_unpaired_entries.begin();

    while (iter_next != this->yet_unpaired_entries.end())
    {

        unpaired_iter = iter_next;
        SamLine * cur_read = const_cast<SamLine *>(*unpaired_iter);
        ++iter_next;

        bool flip_query_strand_flag = this->output_pairs_as_same_strand &&
            cur_read->second_read_in_pair();
        
        if (ignore_bound || less_sam(cur_read, this->low_bound))
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
}
