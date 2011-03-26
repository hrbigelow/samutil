#include "simulation.h"
#include "nucleotide_stats.h"
#include "tools.h"
#include "stats_tools.h"

#include <cmath>
#include <cstdlib>


//parse a file representing the set of possible locus base compositions,
//together with their relative expected proportions.
//line format:  fractionA fractionC fractionG fractionT relative_count
//any number of input lines acceptable
double * parse_basecomp_prior_file(char const* basecomp_prior_file,
                                 double ** base_comps,
                                 double ** base_comp_dist)
{

    size_t num_entries;
    double * buf = ParseNumbersFile(basecomp_prior_file, & num_entries);

    if (num_entries % 5 != 0)
    {
        fprintf(stderr, 
                "Should be 5 entries per line with the format:\n"
                "fractionA fractionC fractionG fractionT locus_fraction_or_count\n"
                "all numbers may be integer or float\n");
        exit(1);
    }
    size_t num_lines = num_entries / 5;

    double * buf2 = new double[num_entries];
    *base_comps = buf2;
    *base_comp_dist = buf2 + (num_lines * 4);

    for (size_t p = 0; p != num_lines; ++p)
    {
        double * src = buf + (p * 5);
        std::copy(src, src + 4, *base_comps + (p * 4));
        std::copy(src + 4, src + 5, *base_comp_dist + p);
    }

    normalize(*base_comp_dist, num_lines, *base_comp_dist);

    delete buf;
    return buf2;

}


//random sample from a stats object to a specified depth, and from
//a specified composition, creating a new nucleotide stats object
LocusSummary sample_locus_from_stats(gsl_rng * rand_gen,
                                     NucleotideStats const& source_stats,
                                     double const* locus_base_comp,
                                     size_t depth)
{
    LocusSummary sampled_locus(source_stats.num_distinct_data, 
                               "", 0, 'A', depth,
                               source_stats.index_mapping);

    //since we don't know ahead of time how many different instances of
    //data we will sample, there is no opportunity for compressing them
    //so, initialize stats index to be 1-to-1 with raw index
    for (size_t ri = 0; ri != sampled_locus.num_distinct_data; ++ri)
    {
        sampled_locus.stats_index[ri] = ri;
    }

    for (size_t d = 0; d != depth; ++d)
    {

        size_t fbase_index =
            SampleDiscreteDistribution(rand_gen, locus_base_comp, 1.0);

        size_t datum_index =
            SampleDiscreteDistribution
            (rand_gen, 
             source_stats.founder_base_likelihood[fbase_index],
             1.0);
        size_t raw_index = datum_index; // because it's 1-to-1 mapping
        
        sampled_locus.raw_counts[raw_index]++;
    }
    
    return sampled_locus;
}
