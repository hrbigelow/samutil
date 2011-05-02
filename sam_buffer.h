#ifndef _SAM_BUFFER_H
#define _SAM_BUFFER_H

#include <set>
#include <utility>

#include "sam_helper.h"
/*
  class for accepting input of SAMLines and outputting them sorted and
  unique.

  a pair of read pair alignments is distinguished by its ordering, based on:

  [contig1, pos1, strand1, cigar1, contig2, pos2, strand2, cigar2]

  if a pair of read pairs have these identical 8 fields, they are
  considered duplicates.  So, when and if a second pair is attempted
  to be inserted into the buffer, the insert will fail.

  
 */



//sort on [rname, pos, strand, cigar, (optionally) qname]
struct LessSAMLinePtr
{
    SAM_ORDER sam_order;
    LessSAMLinePtr(SAM_ORDER _so = SAM_POSITION_RID) : sam_order(_so) { }
    bool operator()(SamLine const* a, SamLine const* b);
};


//sort on [qname]
struct LessSAMLineFragmentIDPtr
{
    bool operator()(SamLine const* a, SamLine const* b);
};


typedef std::pair<SamLine const*, SamLine const*> SAMPTR_PAIR;

//sort on [rname1, pos1, strand1, cigar1, (and optionally) qname1, 
//         rname2, pos2, strand2, cigar2, (and optionally) qname2]
struct LessSAMLinePair
{
    SAM_ORDER sam_order;
    //if true, read pairs with identical alignment coordinates,
    //but different read ids, are considered different.

    LessSAMLinePair(SAM_ORDER _so = SAM_POSITION_RID) : sam_order(_so) { }
    bool operator()(SAMPTR_PAIR const& a, SAMPTR_PAIR const& b);
};


typedef std::set<SAMPTR_PAIR, LessSAMLinePair> PAIRED_READ_SET;
typedef std::pair<PAIRED_READ_SET::iterator, bool> PAIRED_READ_INS;

typedef std::set<SamLine const*, LessSAMLinePtr> SINGLE_READ_SET;
typedef std::pair<SINGLE_READ_SET::iterator, bool> SINGLE_READ_INS;

typedef std::multiset<SamLine const*, LessSAMLineFragmentIDPtr> SINGLE_READ_ID_MSET;


class SamBuffer
{
    SAM_ORDER sam_order;
    bool output_pairs_as_same_strand;
    bool output_is_ones_based;
    bool ignore_duplicate_mapped_pairs;

    //if known, this marks the lowest possible contig position that any 
    //new read will have
    LessSAMLinePtr less_ordering;

 public:

    SamLine const* low_bound;

    PAIRED_READ_SET unique_entry_pairs;
    /* SINGLE_READ_SET single_entries_from_pairs; */

    //entries that are claimed to have a mapped mate, but we haven't seen it yet
    SINGLE_READ_ID_MSET yet_unpaired_entries;

    //entries whose mate is unmapped
    SINGLE_READ_SET unpaired_entries;

    SamBuffer(SAM_ORDER sam_order,
              bool output_pairs_as_same_strand,
              bool output_is_ones_based,
              bool ignore_duplicate_mapped_pairs);

    bool safe_advance_lowbound(SamLine const* _proposed_new_low_bound);

    void update_lowbound(SamLine const* _low_bound);

    //insert an entry, checking whether it is a duplicate.
    bool insert(SamLine const* entry);

    //print entries preceding low_bound, and delete them from memory
    //if 'ignore_bound' set, purge all remaining entries
    void purge(FILE * output_sam_fh, 
               FILE * output_first_fastq_fh, 
               FILE * output_second_fastq_fh, 
               bool ignore_bound);

};

#endif // _SAM_BUFFER_H
