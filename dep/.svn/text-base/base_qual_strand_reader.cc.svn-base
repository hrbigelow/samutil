#include "base_qual_strand_reader.h"
#include "nucleotide_stats.h"
#include "tools.h"
#include "pileup_tools.h"
#include "error_estimate.h"
#include "stats_tools.h"

#include <cstdlib>
#include <cstring>
#include <set>
#include <algorithm>
#include <numeric>


bool BaseQualStrandReader::Datum::operator<(BaseQualStrandReader::Datum const& b) const
{
    return 
        this->called_base < b.called_base 
        || (this->called_base == b.called_base
            && (this->quality < b.quality
                || (this->quality == b.quality
                    && (this->strand < b.strand))));
}


std::string BaseQualStrandReader::Datum::name() const
{
    char dname[50];
    sprintf(dname, "%c_%zi_%c", this->called_base, this->quality, this->strand);
    return std::string(dname);
}


BaseQualStrandReader::Datum 
BaseQualStrandReader::get_datum_from_name(std::string const& name) const
{
    BaseQualStrandReader::Datum d;
    sscanf(name.c_str(), "%c_%zi_%c", &d.called_base, &d.quality, &d.strand);
    return d;
}


//load (base, qual, strand) (ObservedData) tuple data from a locus file,
// checking it against the prior_jpd data for validity
// !!! somehow implement a quality filtering that doesn't break polymorphism...
// ignores N basecalls
std::map<std::string, size_t> parse_locus(PileupSummary const& pileup,
                                          size_t min_quality_score)
{
    
    //create a temporary structure of raw counts using 4 x Q x 2 flat array dimension,
    //Q is highest quality + 1
    std::map<std::string, size_t>::const_iterator name_iter;
    size_t highest_quality = 0;
    for (size_t qi = 0; qi != pileup._read_depth; ++qi)
    {
        highest_quality = 
            std::max(highest_quality, 
                     QualityCodeToQuality(pileup._quality_codes[qi]));
    }

    size_t const num_s = 2;
    size_t const num_qs = num_s * (highest_quality + 1);
    size_t const num_bqs = 4 * num_qs;
    size_t * raw_counts_flat = new size_t[num_bqs];
    std::fill(raw_counts_flat, raw_counts_flat + num_bqs, 0);

    size_t reads_counted = 0;

    for (size_t r = 0; r < static_cast<size_t>(pileup._read_depth); ++r)
    {
        size_t quality = QualityCodeToQuality(pileup._quality_codes[r]);
        char basecall = pileup._bases[r];
        size_t basecall_index = Nucleotide::base_to_index[static_cast<size_t>(basecall)];
        if (basecall_index >= 4)
        {
            //since 'N' is common in pileup, but meaningless, we ignore it silently
            continue;
        }
        if (quality < min_quality_score)
        {
            continue;
        }
        char strand = isupper(basecall) ? '+' : '-';
        size_t strand_index = strand == '+' ? 0 : 1;
        size_t flat_index = (basecall_index * num_qs) + (quality * num_s) + strand_index;
        raw_counts_flat[flat_index]++;
        ++reads_counted;
    }

    std::map<std::string, size_t> counts;

    for (size_t flat_index = 0; flat_index != num_bqs; ++flat_index)
    {
        if (raw_counts_flat[flat_index] > 0)
        {
            char basecall = Nucleotide::bases_upper[flat_index / num_qs];
            size_t quality = (flat_index % num_qs) / num_s; 
            char strand = (flat_index % num_s) == 0 ? '+' : '-';

            BaseQualStrandReader::Datum datum = 
                { 
                    basecall, quality, strand, 0, 0
                };
        
            if (counts.find(datum.name()) != counts.end())
            {
                fprintf(stderr, "BaseQualStrandData::parse_locus: unknown error.\n");
                exit(13);
            }
            counts[datum.name()] = raw_counts_flat[flat_index];
        }
    }
    delete raw_counts_flat;
    return counts;

}



LocusSummary 
BaseQualStrandReader::get_next_locus(NucleotideStats const& nuc_stats,
                                     void const* extra)
{

    size_t min_quality_score = * static_cast<size_t const*>(extra);

    PileupSummary pileup(0);
    bool succeeded = pileup.load_line(this->fh);
    
    if (! succeeded)
    {
        fprintf(stderr, "Couldn't parse pileup line\n");
        exit(1);
    }
    
    std::map<std::string, size_t> counts = 
        parse_locus(pileup, min_quality_score);


    std::map<std::string, size_t>::const_iterator prior_iter;

    std::map<std::string, size_t>::const_iterator cit;
    size_t reads_counted = 0;
    for (cit = counts.begin(); cit != counts.end(); ++cit)
    {
        reads_counted += (*cit).second;
    }

    LocusSummary locus(counts.size(), pileup._reference, pileup._position,
                       pileup._reference_base, reads_counted,
                       nuc_stats.index_mapping);
    
    //check for presence in name_map
    size_t raw_index;
    for (cit = counts.begin(), raw_index = 0; 
         cit != counts.end(); 
         ++cit, ++raw_index)
    {
        prior_iter = nuc_stats.name_mapping.find((*cit).first);
        if (prior_iter == nuc_stats.name_mapping.end())
        {
            //this should never happen
            fprintf(stderr, "BaseQualStrandReader::get_next_locus: "
                    "found datum %s not appearing in prior data.\n",
                    (*cit).first.c_str());

            exit(10);
        }
        size_t stat_index = (*prior_iter).second;
        locus.raw_counts[raw_index] = (*cit).second;
        locus.stats_index[raw_index] = stat_index;
    }    

    return locus;

}


void
BaseQualStrandReader::compute_strand_marginal(JPD_DATA const& counts_map,
                                              double * pos_strand_counts,
                                              double * neg_strand_counts) const
{
    *pos_strand_counts = 0.0;
    *neg_strand_counts = 0.0;

    JPD_DATA::const_iterator cit;
    for (cit = counts_map.begin(); cit != counts_map.end(); ++cit)
    {
        Datum datum = this->get_datum_from_name((*cit).first);
        double & counts = 
            datum.strand == '+' ? *pos_strand_counts : *neg_strand_counts;

        counts += (*cit).second.total();
    }
}


//return a new JPD such that the strand marginals become those specified
//total mass remain the same as original counts
JPD_DATA
BaseQualStrandReader::normalize_strand_marginal(double const* desired_strand_marginal,
                                                JPD_DATA const& counts_map) const
{
    

    double current_strand_marginal[2];
    this->compute_strand_marginal(counts_map, current_strand_marginal,
                                  current_strand_marginal + 1);

    normalize(current_strand_marginal, 2, current_strand_marginal);
    double sum_of_desired = 
        desired_strand_marginal[0]
        + desired_strand_marginal[1];

    double flow_ratio_marginal[2];
    flow_ratio_marginal[0] = 
        desired_strand_marginal[0] 
        / (current_strand_marginal[0] * sum_of_desired);

    flow_ratio_marginal[1] = 
        desired_strand_marginal[1] 
        / (current_strand_marginal[1] * sum_of_desired);

    assert(! isnan(flow_ratio_marginal[0]));
    assert(! isnan(flow_ratio_marginal[1]));

    JPD_DATA new_counts_map(counts_map);

    JPD_DATA::iterator cit;
    for (cit = new_counts_map.begin(); cit != new_counts_map.end(); ++cit)
    {
        Datum datum = this->get_datum_from_name((*cit).first);
        double & flow = 
            datum.strand == '+' ? flow_ratio_marginal[0] : flow_ratio_marginal[1];

        assert(! isnan(flow));

        (*cit).second.multiply(flow);
    }

//     double new_strand_marginal[2];
//     this->compute_strand_marginal(new_counts_map, new_strand_marginal,
//                                   new_strand_marginal + 1);

//     assert(new_strand_marginal[0] == desired_strand_marginal[0]);
//     assert(new_strand_marginal[1] == desired_strand_marginal[1]);

    return new_counts_map;
}


//make a LocusSummary identical to the one from full_locus, but zero
//out all counts for the other strand.
LocusSummary 
BaseQualStrandReader::locus_data_by_strand(LocusSummary const& full_locus, 
                                           char strand) const
{
    //make an identical copy
    LocusSummary strand_locus(full_locus);

    for (size_t raw_index = 0; raw_index != strand_locus.num_distinct_data; 
         ++raw_index)
    {
        size_t datum_index = strand_locus.stats_index[raw_index];
        BaseQualStrandReader::Datum datum = 
            get_datum_from_name(strand_locus.index_mapping[datum_index]);
        if (datum.strand != strand)
        {
            strand_locus.read_depth -= strand_locus.raw_counts[raw_index];
            strand_locus.raw_counts[raw_index] = 0.0;
        }
    }
    return strand_locus;
}
