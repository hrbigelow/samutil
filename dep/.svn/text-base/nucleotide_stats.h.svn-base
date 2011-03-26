#ifndef _NUCLEOTIDE_STATS_H
#define _NUCLEOTIDE_STATS_H

#include <map>
#include <string>

namespace Nucleotide 
{
    //transforms ACGT and acgt to 0123, everything else to 4
    extern int const base_to_index[];
    extern char const* bases_upper;
};


struct nuc_frequency
{
    double data[4];
    nuc_frequency(double const _data[4])
    {
        std::copy(_data, _data + 4, this->data);
    }

    nuc_frequency(nuc_frequency const& a)
    {
        std::copy(a.data, a.data + 4, this->data);
    }

    double total() const
    {
        return data[0] + data[1] + data[2] + data[3];
    }
    void multiply(double factor)
    {
        data[0] *= factor;
        data[1] *= factor;
        data[2] *= factor;
        data[3] *= factor;
    }
};
    

typedef std::map<std::string, nuc_frequency> JPD_DATA;

class LocusSummary;


/* 
 */
class NucleotideStats {

    double * jpd_buffer;
    double * cpd_buffer;

 public:
    double founder_base_marginal[4];
    double * complete_jpd[4];
    double * founder_base_likelihood[4];
    size_t num_distinct_data;
    std::map<std::string, size_t> name_mapping;
    std::string * index_mapping;

    NucleotideStats(size_t num_distinct_data_types);
    ~NucleotideStats();
    void initialize(JPD_DATA const& counts_map);
    JPD_DATA make_per_locus_stats(LocusSummary const& locus);

};


class LocusSummary {

 public:
    LocusSummary(size_t _num_distinct_data,
                 char const* _reference, size_t _position,
                 char _reference_base, size_t _read_depth,
                 std::string const* _index_mapping);

    LocusSummary(LocusSummary const&);
    ~LocusSummary();

    double * raw_counts;

    //maps index in LocusSummary::raw_counts to index in
    //NucleotideStats::founder_base_likelihood
    size_t * stats_index;

    size_t num_distinct_data;
    char reference[100];
    size_t position;
    char reference_base;
    size_t read_depth;

    //points to the index_mapping in NucleotideStats.  not owned here
    std::string const* index_mapping;
};

#endif // _NUCLEOTIDE_STATS_H
