#ifndef _BASE_QUAL_STRAND_READER_H
#define _BASE_QUAL_STRAND_READER_H

#include "nucleotide_reader.h"


class BaseQualStrandReader : public NucleotideReader
{
 public:

    struct Datum {
        char called_base;
        size_t quality;
        char strand;
        size_t order_index;
        double frequency;
        std::string name() const;

        bool operator<(Datum const& b) const;
    };


    BaseQualStrandReader() : NucleotideReader() { }
        LocusSummary get_next_locus(NucleotideStats const& nuc_stats, void const* extra);

    Datum get_datum_from_name(std::string const& name) const;
    void compute_strand_marginal(JPD_DATA const& counts_map, 
                                 double *pos_strand_counts,
                                 double *neg_strand_counts) const;

    JPD_DATA normalize_strand_marginal(double const* desired_strand_marginal,
                                       JPD_DATA const& counts_map) const;

    LocusSummary locus_data_by_strand(LocusSummary const& full_locus, char strand) const;

};


#endif // _BASE_QUAL_STRAND_READER_H
