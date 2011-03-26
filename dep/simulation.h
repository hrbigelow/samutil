#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <utility>

#include <gsl/gsl_rng.h>

class ErrorEstimate;
class LocusSummary;
class NucleotideStats;

//sets the output distribution to uniform for all non-major indices,
//and m for the i.e {0.6, 0.1, 0.1, 0.1, 0.1}
void SetMajorMinorComposition(double major_prob, int major_index,
                              double * output_distribution, int size);

std::pair<char, int> SimulateBaseMeasure(gsl_rng * rand_gen,
                                         int const* quality_counts, 
                                         int total_counts, char base);

double * parse_basecomp_prior_file(char const* basecomp_prior_file,
                                   double ** base_comps,
                                   double ** base_comp_dist);

LocusSummary sample_locus_from_stats(gsl_rng * rand_gen,
                                     NucleotideStats const& source_stats,
                                     double const* locus_base_comp,
                                     size_t depth);

#endif // _SIMULATION_H
