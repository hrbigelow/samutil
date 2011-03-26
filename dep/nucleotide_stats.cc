#include "nucleotide_stats.h"
#include "stats_tools.h"

#include <cstring>
#include <cassert>
#include <numeric>
#include <cmath>


NucleotideStats::NucleotideStats(size_t _num_distinct_data) :
    num_distinct_data(_num_distinct_data)
{
    size_t D = _num_distinct_data;

    this->jpd_buffer = new double[D * 4];
    this->cpd_buffer = new double[D * 4];

    for (size_t b = 0; b != 4; ++b)
    {
        this->complete_jpd[b] = this->jpd_buffer + (b * D);
        this->founder_base_likelihood[b] = this->cpd_buffer + (b * D);
    }
    this->index_mapping = new std::string[D];
}

NucleotideStats::~NucleotideStats()
{
    delete this->jpd_buffer;
    delete this->cpd_buffer;
    delete[] this->index_mapping;
}

namespace Nucleotide
{
    int const base_to_index[] =
        {
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
            4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
            4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
            4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
        };
    char const* bases_upper = "ACGTN";
};


//initialize all distributions from the counts_map.  counts_map
//is not necessarily normalized.  any entries in this stats object
//without corresponding entries in counts_map are regarded as having
//zero counts
void NucleotideStats::initialize(JPD_DATA const& counts_map)
{
    size_t D = this->num_distinct_data;

    JPD_DATA::const_iterator cit;

    std::fill(this->jpd_buffer, this->jpd_buffer + (4 * D), 0.0);

    size_t datum_index;
    for (cit = counts_map.begin(), datum_index = 0; 
         cit != counts_map.end(); 
         ++cit, ++datum_index)
    {
        this->index_mapping[datum_index] = (*cit).first;
        this->name_mapping[(*cit).first] = datum_index;
        for (size_t b = 0; b != 4; ++b)
        {
            this->complete_jpd[b][datum_index] = (*cit).second.data[b];
        }
    }
    
    normalize(this->jpd_buffer, 4 * D, this->jpd_buffer);
    
    for (size_t b = 0; b != 4; ++b)
    {
        this->founder_base_marginal[b] =
            std::accumulate(this->complete_jpd[b],
                            this->complete_jpd[b] + D, 0.0);
    }

    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t di = 0; di != D; ++di)
        {
            this->founder_base_likelihood[b][di] =
                this->complete_jpd[b][di]
                / this->founder_base_marginal[b];
        }
    }

}


//creates a new NucleotideStats object with the same data set
//as the calling one, but with a marginal distribution of observed data
//reflected by 'locus', but the shape P(fb|obs) of the calling object
JPD_DATA
NucleotideStats::make_per_locus_stats(LocusSummary const& locus)
{
    JPD_DATA counts_map;
    for (size_t di = 0; di != this->num_distinct_data; ++di)
    {
        double zero_counts[] = { 0, 0, 0, 0 };
        counts_map.insert(std::make_pair(this->index_mapping[di], nuc_frequency(zero_counts)));
    }

    double locus_slice[4];
    for (size_t raw_index = 0; raw_index != locus.num_distinct_data; ++raw_index)
    {
        size_t di = locus.stats_index[raw_index];
        double orig_marginal = 
            this->complete_jpd[0][di]
            + this->complete_jpd[1][di]
            + this->complete_jpd[2][di]
            + this->complete_jpd[3][di];

        double adjust_factor;
        if (orig_marginal == 0)
        {
            assert(locus.raw_counts[raw_index] == 0);
            adjust_factor = 1.0;
        }
        else
        {
            adjust_factor = locus.raw_counts[raw_index] / orig_marginal;
        }

        //double adjust_factor = 1.0;
        for (size_t bi = 0; bi != 4; ++bi)
        {
            locus_slice[bi] = this->complete_jpd[bi][di] * adjust_factor;
        }
        counts_map.erase(this->index_mapping[di]);
        counts_map.insert(std::make_pair(this->index_mapping[di], nuc_frequency(locus_slice)));
    }
    return counts_map;
}


LocusSummary::LocusSummary(size_t _num_distinct_data,
                           char const* _reference, size_t _position,
                           char _reference_base, size_t _read_depth,
                           std::string const* _index_mapping) : 
    num_distinct_data(_num_distinct_data),
    position(_position), 
    reference_base(_reference_base), 
    read_depth(_read_depth),
    index_mapping(_index_mapping)
{
    strcpy(reference, _reference);
    raw_counts = new double[_num_distinct_data];
    stats_index = new size_t[_num_distinct_data];

    std::fill(raw_counts, raw_counts + _num_distinct_data, 0.0);
}


LocusSummary::LocusSummary(LocusSummary const& ls)
    : num_distinct_data(ls.num_distinct_data),
      position(ls.position),
      reference_base(ls.reference_base),
      read_depth(ls.read_depth),
      index_mapping(ls.index_mapping)
{
    strcpy(this->reference, ls.reference);
    this->raw_counts = new double[this->num_distinct_data];
    this->stats_index = new size_t[this->num_distinct_data];
    std::copy(ls.raw_counts, ls.raw_counts + ls.num_distinct_data, this->raw_counts);
    std::copy(ls.stats_index, ls.stats_index + ls.num_distinct_data, this->stats_index);
}


LocusSummary::~LocusSummary()
{
    delete raw_counts;
    delete stats_index;
}
