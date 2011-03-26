#include "alignment_stats.h"

#include <cassert>
#include <cmath>

RefStrandPair::RefStrandPair(char const* _guide_ref, Strand _guide_strand,
                             char const* _test_ref, Strand _test_strand) :
    guide_strand(_guide_strand), test_strand(_test_strand)
{
    strcpy(guide_ref, _guide_ref);
    strcpy(test_ref, _test_ref);
}


bool 
operator<(RefStrandPair const& a, RefStrandPair const& b)
{
    int guide_cmp = strcmp(a.guide_ref, b.guide_ref);
    int test_cmp = strcmp(a.test_ref, b.test_ref);
    
    return guide_cmp < 0
        || (guide_cmp == 0
            && (a.guide_strand < b.guide_strand
                || (a.guide_strand == b.guide_strand
                    && (test_cmp < 0
                        || (test_cmp == 0
                            && (a.test_strand < b.test_strand))))));
}


bool 
operator<(AlignCompareStats const& a, AlignCompareStats const& b)
{
    return a.num_correct < b.num_correct
        || (a.num_correct == b.num_correct
            && (a.num_overplaced < b.num_overplaced
                || (a.num_overplaced == b.num_overplaced
                    && (a.num_underplaced < b.num_underplaced
                        || (a.num_underplaced == b.num_underplaced
                            && (a.num_misplaced < b.num_misplaced))))));
}


//output the bin number from a coordinate
int OffsetHistogram::bin(int64_t coordinate) const
{
    return 
        static_cast<int>(floor(static_cast<float>(coordinate) / 
                               static_cast<float>(this->binsize)));
}


int64_t OffsetHistogram::lowbin_bound(int bin) const
{
    return static_cast<int64_t>(bin) * this->binsize;
}


int64_t OffsetHistogram::highbin_bound(int bin) const
{
    return static_cast<int64_t>(bin + 1) * this->binsize - 1;
}



AlignCompareStats 
TallyAlignmentStats(SamLine const& guide,
                    SamLine const& test,
                    std::vector<Cigar::UnitComparison> const& cigar_comparison,
                    size_t histo_binsize,
                    std::map<RefStrandPair, OffsetHistogram> * misplaced_bases)
{
    AlignCompareStats align_compare_stats(0,0,0,0);

    for (size_t c = 0; c != cigar_comparison.size(); ++c)
    {
        Cigar::UnitComparison const& cc = cigar_comparison[c];
                
        bool placed_on_guide, placed_on_test;

        switch (cc.guide_op)
        {
        case Cigar::M:
        case Cigar::I:
            placed_on_guide = true;
            break;
        case Cigar::D:
        case Cigar::N:
        case Cigar::S:
        case Cigar::H:
        case Cigar::P:
            placed_on_guide = false;
            break;

        case Cigar::None:
            assert(false);
            break;
        }

        switch (cc.test_op)
        {
        case Cigar::M:
        case Cigar::I:
            placed_on_test = true;
            break;
        case Cigar::D:
        case Cigar::N:
        case Cigar::S:
        case Cigar::H:
        case Cigar::P:
            placed_on_test = false;
            break;

        case Cigar::None:
            assert(false);
            break;
        }

        if (placed_on_guide)
        {
            if (placed_on_test)
            {
                Strand guide_strand = guide.query_on_pos_strand() ? POS : NEG;
                Strand test_strand = test.query_on_pos_strand() ? POS : NEG;

                if (cc.offset == 0
                    && (strcmp(guide.rname, test.rname) == 0)
                    && (guide_strand == test_strand)
                    )
                {
                    //correctly placed
                    align_compare_stats.num_correct += cc.length;
                }
                else
                {
                    std::map<RefStrandPair, OffsetHistogram> & mbc = (*misplaced_bases);
                    RefStrandPair rp(guide.rname, guide_strand,
                                     test.rname, test_strand);
                    
                    if (mbc.find(rp) == mbc.end())
                    {
                        mbc.insert(std::make_pair(rp, OffsetHistogram(histo_binsize)));
                    }
                    
                    if (mbc[rp].counts.find(cc.offset) == mbc[rp].counts.end())
                    {
                        mbc[rp].counts.insert(std::make_pair(cc.offset, 0));
                    }
                    mbc[rp].counts[cc.offset] += cc.length;

                    align_compare_stats.num_misplaced += cc.length;
                }
            }
            else
            {
                align_compare_stats.num_underplaced += cc.length;
            }
        }
        else
        {
            //not placed on guide
            if (placed_on_test)
            {
                align_compare_stats.num_overplaced += cc.length;
            }
        }
    }
    return align_compare_stats;
}


//compare the set of cigar operations and offsets, loading them into a
//comparison object
std::vector<Cigar::UnitComparison> 
TallyCigarOffsets(Cigar::CIGAR_VEC guide_cigar,
                  size_t guide_pos,
                  Cigar::CIGAR_VEC test_cigar,
                  size_t test_pos)
{
    size_t gpos = guide_pos;
    size_t tpos = test_pos;
    size_t gchunk = 0;
    size_t tchunk = 0;

    size_t step;

    std::vector<Cigar::UnitComparison> comparison;

    //slight overkill, but might be about this size
    comparison.reserve(guide_cigar.size() + test_cigar.size());

    while (gchunk < guide_cigar.size()
           && tchunk < test_cigar.size())
    {

        if (guide_cigar[gchunk].op == Cigar::D ||
            guide_cigar[gchunk].op == Cigar::N)
        {
            gpos += guide_cigar[gchunk].length;
            ++gchunk;
        }

        if (test_cigar[tchunk].op == Cigar::D ||
            test_cigar[tchunk].op == Cigar::N)
        {
            tpos += test_cigar[tchunk].length;
            ++tchunk;
        }

        Cigar::Unit & guide_unit = guide_cigar[gchunk];
        Cigar::Unit & test_unit = test_cigar[tchunk];

        step = std::min(guide_unit.length, test_unit.length);

        Cigar::OpComp op_pair = static_cast<Cigar::OpComp>(guide_unit.op * 7 + test_unit.op);

        Cigar::UnitComparison cigar_comp = { op_pair, guide_unit.op, test_unit.op, step, gpos - tpos };
        comparison.push_back(cigar_comp);

        guide_unit.length -= step;
        test_unit.length -= step;

        gpos += step;
        tpos += step;

        if (guide_unit.length == 0)
        {
            ++gchunk;
        }
        if (test_unit.length == 0)
        {
            ++tchunk;
        }
    }
    assert(gchunk == guide_cigar.size());
    assert(tchunk == test_cigar.size());

    return comparison;
}


char Bool2Strand(bool pos_if_true)
{
    return pos_if_true ? 'P' : 'N';
}
