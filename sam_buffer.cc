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
