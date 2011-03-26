#ifndef _NUCLEOTIDE_READER_H
#define _NUCLEOTIDE_READER_H

#include "nucleotide_stats.h"

#include <cstdio>
#include <map>
#include <string>

/* class LocusSummary; */
class NucleotideStats;

class NucleotideReader {

protected:
    FILE * fh;

public:
    
    NucleotideReader();
    ~NucleotideReader();
    bool initialize(char const* file);
    NucleotideStats read_from_rdb(char const* rdb_file);

    virtual LocusSummary 
        get_next_locus(NucleotideStats const& nuc_stats, void const* extra) = 0;
    
    virtual void
        compute_strand_marginal(JPD_DATA const& counts_map, 
                                double *pos_strand_counts,
                                double *neg_strand_counts) const = 0;
    
    virtual JPD_DATA
        normalize_strand_marginal(double const* desired_strand_marginal,
                                  JPD_DATA const& counts_map) const = 0;

    virtual LocusSummary
        locus_data_by_strand(LocusSummary const& full_locus, char strand) const = 0;
    
    //tally the total frequency of positive and negative stranded reads
    bool more_loci();
};


JPD_DATA parse_jpd_rdb_file(char const* rdb_file);
    

#endif // _NUCLEOTIDE_READER_H
