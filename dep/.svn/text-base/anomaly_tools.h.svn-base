#ifndef _ANOMALY_TOOLS_H
#define _ANOMALY_TOOLS_H

class Posterior;
class NucleotideReader;
class LocusSummary;

#include "nucleotide_stats.h"

#include <cstddef>

double williams_moment_match_ratio(NucleotideStats const& stats,
                                   size_t N,
                                   double const* evaluation_point);


double strand_locus_anomaly_score(Posterior & posterior,
                                  JPD_DATA const& global_counts,
                                  LocusSummary const& full_locus,
                                  NucleotideReader const* data_reader,
                                  char strand,
                                  bool verbose);


#endif // _ANOMALY_TOOLS_H
