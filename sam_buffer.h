#ifndef _SAM_BUFFER_H
#define _SAM_BUFFER_H

#include <set>
#include <utility>

#include "sam_helper.h"
#include "sam_order.h"

/*
The SAM buffer makes up for a shortcoming in the SAM spec that the two
(or in the future more than two) reads generated from a single
fragment of DNA are stored as separate records, when they should
instead be stored as a single one.

So, the SAM buffer is 'fragment-centric', where it's concept of a
'fragment' is currently one of:

1. the pair of SAM records that are paired-in-sequencing and
mapped-in-proper-pair

2. individual SAM records that are not paired-in-sequencing.

The SAM buffer assumes that the alignments it receives are 'valid' in
the sense that if a fragment is paired-in-sequencing, then it is also
either entirely unmapped, or mapped-in-proper-pair.

What it does:

1. collects individual SAM records in an arbitrary order, and allows
outputting them in a specified order between pairs, but keeping paired
records together.

2. de-duplicates records or record-pairs, based on specified ordering
criterion.

3. checks the pair-validity of the SAM file, that records claiming to
have a pair actually do.

The recommended practice is to use the SAM buffer to load, pair and
convert ordinary SAM records to rSAM. To avoid overflowing memory, the
buffer must be emptied periodically. However, it can only be emptied
when all entries are properly paired.

*/


//sort on [rname, pos, strand, cigar, (optionally) qname]
struct LessSAMLinePtr
{
    SamOrder const* sam_order;
    LessSAMLinePtr(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b) const;
};


//sort on [qname]
struct LessSAMLineFragmentIDPtr
{
    SamOrder const* sam_order;
    LessSAMLineFragmentIDPtr(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b) const;
};


struct LessSAMLinePtrMatePair
{
    SamOrder const* sam_order;
    LessSAMLinePtrMatePair(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b) const;
};



typedef std::set<SamLine const*, LessSAMLinePtr> SINGLE_READ_SET;
typedef std::pair<SINGLE_READ_SET::iterator, bool> SINGLE_READ_INS;

typedef std::multiset<SamLine const*, LessSAMLineFragmentIDPtr> SINGLE_READ_ID_MSET;


struct InsertResult
{
    SINGLE_READ_SET::iterator surviving_entry;
    bool was_inserted;
    SamLine * remaining_entry;
    InsertResult(SINGLE_READ_SET::iterator _iter, bool _wi, SamLine * _sam);
    InsertResult();
};


class SamBuffer
{
    bool output_pairs_as_same_strand;

    //if known, this marks the lowest possible contig position that any 
    //new read will have
    LessSAMLinePtr less_ordering;
    LessSAMLinePtrMatePair less_entry_matepair;

 public:

    SamOrder const* sam_order;

    SINGLE_READ_SET unique_entries;
    /* SINGLE_READ_SET single_entries_from_pairs; */

    //entries that are claimed to have a mapped mate, but we haven't seen it yet
    SINGLE_READ_ID_MSET incomplete_entries;

    SamBuffer(SamOrder const* sam_order, bool output_pairs_as_same_strand);

    // insert an entry, checking whether it is a duplicate.
    // this follows the same logic as std::set::insert() in its return type.
    InsertResult insert(SamLine * entry);

    void replace(SINGLE_READ_SET::iterator iter,
                 SamLine * replacement_entry);

};

#endif // _SAM_BUFFER_H
