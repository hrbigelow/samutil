#ifndef _ALIGNMENT_STATS_H
#define _ALIGNMENT_STATS_H

#include "sam_helper.h"
#include "cigar_ops.h"

#include <vector>
#include <utility>
#include <map>

struct RefStrandPair
{
    char guide_ref[100];
    Strand guide_strand;
    char test_ref[100];
    Strand test_strand;

    RefStrandPair(char const* guide_ref, Strand guide_strand,
                  char const* test_ref, Strand test_strand);

    bool friend operator<(RefStrandPair const& a, RefStrandPair const& b);

};

/* struct misplaced_block */
/* { */
/*     char guide_ref[100]; */
/*     Strand guide_strand; */
/*     char test_ref[100]; */
/*     Strand test_strand; */
/*     int64_t offset; */
/*     size_t block_length; */
    
/*     misplaced_block(char * guide_ref, Strand guide_strand, */
/*                     char * test_ref, Strand test_strand,  */
/*                     int64_t o, size_t l); */
/* }; */


//summary statistics describing accuracy of an alignment
struct AlignCompareStats
{
    size_t num_correct;
    size_t num_overplaced;
    size_t num_underplaced;
    size_t num_misplaced;
    size_t num_total() const
    {
        return num_correct + num_overplaced + num_underplaced + num_misplaced;
    }

    AlignCompareStats(size_t c, size_t o, size_t u, size_t m) :
        num_correct(c), num_overplaced(o), num_underplaced(u),
        num_misplaced(m) { }

    bool friend operator<(AlignCompareStats const& a, AlignCompareStats const& b);
};


typedef std::map<int64_t, size_t> OFFSET_COUNTS;

struct OffsetHistogram
{
    size_t binsize;
    OFFSET_COUNTS counts;

    OffsetHistogram(size_t b) : binsize(b) { }
    OffsetHistogram() : binsize(0) { }

    int bin(int64_t coordinate) const;
    int64_t lowbin_bound(int bin) const;
    int64_t highbin_bound(int bin) const;
};


AlignCompareStats 
TallyAlignmentStats(SamLine const& guide,
                    SamLine const& test,
                    std::vector<Cigar::UnitComparison> const& cigar_comparison,
                    size_t histo_binsize,
                    std::map<RefStrandPair, OffsetHistogram> * misplaced_bases);

std::vector<Cigar::UnitComparison> 
TallyCigarOffsets(Cigar::CIGAR_VEC guide_cigar,
                  size_t guide_pos,
                  Cigar::CIGAR_VEC test_cigar,
                  size_t test_pos);


char Bool2Strand(bool pos_if_true);




#endif // _ALIGNMENT_STATS_H
