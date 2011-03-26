#include "pretty_plot.h"

#include "dep/tools.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <limits.h>

size_t pseudo_intron_length_def = 25;


int pretty_plot_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "pretty_plot\n"
            "Author: Henry Bigelow (hbigelow@amgen.com)\n\n"
            "pretty_plot gtf_sam        Condense introns in gtf and sam file\n"
            "pretty_plot graph          Condense and annotate exons and loci for graphing\n\n"
            );
    return 1;
}


int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        return pretty_plot_usage();
    }
    else if (strcmp(argv[1], "gtf_sam") == 0)
    {
        return main_pretty_plot_gtf_sam(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "graph") == 0)
    {
        return main_pretty_plot_graph(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
    
}
