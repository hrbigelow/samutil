#include "align_eval_coverage.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <getopt.h>
#include <cstdlib>

#include "dep/tools.h"

int align_eval_coverage_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval coverage [OPTIONS] jumps.txt coverage.txt\n\n"
            "Options:\n\n"
            "-u  FLAG   (unique) assume the workflow came from unique alignments\n"
            "           if so, 'correct' depth should not exceed 'guide' depth and\n"
            "           is an error if it does [false].\n"
            );

    fprintf(stderr,
            "jumps.txt has fields:\n\n"
            "<contig> <position> <guide_jump> <correct_jump> <error_jump>\n"
            "See 'align_eval raw' to produce jumps.txt file from SAM alignment\n\n");
    return 1;
}


int main_align_eval_coverage(int argc, char ** argv)
{
    char c;
    bool assume_unique_alignments = false;

    while ((c = getopt(argc, argv, "u")) >= 0)
    {
        switch(c)
        {
        case 'u': assume_unique_alignments = true; break;
        default: return align_eval_coverage_usage(); break;
        }
    }
    
    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return align_eval_coverage_usage();
    }

    char const* jumps_file = argv[optind];
    char const* coverage_file = argv[optind + 1];

    FILE * jumps_fh = open_or_die(jumps_file, "r", "Input Jumps file");
    FILE * coverage_fh = open_or_die(coverage_file, "w", "Output coverage file");

    //maintain a running current total of guide, correct, and error.
    char contig[256] = "";
    size_t bound, prev_bound = 0;
    int guide_jump, correct_jump, error_jump;
    int guide = 0, correct = 0, error = 0;

    while (! feof(jumps_fh))
    {
        fscanf(jumps_fh, "%s\t%zu\t%i\t%i\t%i\n",
               contig, &bound, &guide_jump, &correct_jump, &error_jump);

        if (guide != 0 || correct != 0 || error != 0)
        {
            //only print this if one or more values are nonzero
            for (size_t pos = prev_bound; pos != bound; ++pos)
            {
                fprintf(coverage_fh, "%s\t%Zu\t%i\t%i\t%i\n", contig, pos + 1, guide, correct, error);
            }
            if (assume_unique_alignments && guide < correct)
            {
                fprintf(stderr, "Error: guide depth %i less than correct depth %i "
                        "and assuming unique alignments.\n", guide, correct);
                exit(1);
            }
        }

        guide += guide_jump;
        correct += correct_jump;
        error += error_jump;

        prev_bound = bound;
    }

    assert(guide == 0);
    assert(correct == 0);
    assert(error == 0);

    fclose(jumps_fh);
    fclose(coverage_fh);

    return 0;
}
