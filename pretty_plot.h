#ifndef _PRETTY_PLOT_H
#define _PRETTY_PLOT_H

#include "cigar_ops.h"

#include <map>
#include <string>

#include "nclist.h"

int pretty_plot_graph_usage(size_t ldef);
int pretty_plot_gtf_sam_usage(size_t ldef);

int main_pretty_plot_graph(int argc, char **argv);
int main_pretty_plot_gtf_sam(int argc, char **argv);

extern size_t pseudo_intron_length_def;

#endif // _PRETTY_PLOT_H

