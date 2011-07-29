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

For memory efficiency, though, there is a single wrinkle here:

In order not to store the entire set of records processed, the SAM
buffer maintains a marker called a 'low_bound'. The low bound is a
claim from the client that says "I guarantee that I will not insert
any further SAM records that fall below this bound, according to the
specified output ordering.

This claim is useful because it allows the caller to 'purge' the SAM
buffer up to (but not including) this low bound without violating the
output ordering or missing any opportunities to de-duplicate the
records.



  
 */


//sort on [rname, pos, strand, cigar, (optionally) qname]
struct LessSAMLinePtr
{
    SamOrder const* sam_order;
    LessSAMLinePtr(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b);
};


//sort on [qname]
struct LessSAMLineFragmentIDPtr
{
    SamOrder const* sam_order;
    LessSAMLineFragmentIDPtr(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b);
};


typedef std::pair<SamLine const*, SamLine const*> SAMPTR_PAIR;

//sort on [rname1, pos1, strand1, cigar1, (and optionally) qname1, 
//         rname2, pos2, strand2, cigar2, (and optionally) qname2]
struct LessSAMLinePair
{
    SamOrder const* sam_order;
    //if true, read pairs with identical alignment coordinates,
    //but different read ids, are considered different.

    LessSAMLinePair(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SAMPTR_PAIR const& a, SAMPTR_PAIR const& b);
};


struct LessSAMLinePtrMatePair
{
    SamOrder const* sam_order;
    LessSAMLinePtrMatePair(SamOrder const* _so = NULL) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b);
};




typedef std::set<SAMPTR_PAIR, LessSAMLinePair> PAIRED_READ_SET;
typedef std::pair<PAIRED_READ_SET::iterator, bool> PAIRED_READ_INS;

typedef std::set<SamLine const*, LessSAMLinePtr> SINGLE_READ_SET;
typedef std::pair<SINGLE_READ_SET::iterator, bool> SINGLE_READ_INS;

typedef std::multiset<SamLine const*, LessSAMLineFragmentIDPtr> SINGLE_READ_ID_MSET;


class SamBuffer
{
    bool output_pairs_as_same_strand;

    //if known, this marks the lowest possible contig position that any 
    //new read will have
    LessSAMLinePtr less_ordering;
    LessSAMLinePtrMatePair less_entry_matepair;

 public:

    SamOrder const* sam_order;
    SamLine const* low_bound;

    PAIRED_READ_SET unique_entry_pairs;
    /* SINGLE_READ_SET single_entries_from_pairs; */

    //entries that are claimed to have a mapped mate, but we haven't seen it yet
    SINGLE_READ_ID_MSET yet_unpaired_entries;

    //entries whose mate is unmapped
    SINGLE_READ_SET unpaired_entries;

    SamBuffer(SamOrder const* sam_order, bool output_pairs_as_same_strand);

    bool safe_advance_lowbound(SamLine const* _proposed_new_low_bound);

    void update_lowbound(SamLine const* _low_bound);

    //insert an entry, checking whether it is a duplicate.
    std::pair<SamLine const*, bool> insert(SamLine const* entry);

    //print entries preceding low_bound, and delete them from memory
    //if 'ignore_bound' set, purge all remaining entries
    void purge(FILE * output_sam_fh, 
               FILE * output_first_fastq_fh, 
               FILE * output_second_fastq_fh, 
               bool ignore_bound);

};

#endif // _SAM_BUFFER_H
