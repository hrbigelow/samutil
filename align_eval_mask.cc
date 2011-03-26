#include "align_eval_mask.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <getopt.h>
#include <cstdlib>

#include "dep/tools.h"

int align_eval_mask_usage(float cdef, float edef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval mask [OPTIONS] jumps.txt regions.txt\n\n"
            "Options:\n\n"
            "-c  FLOAT   minimum fraction correct bases.  loci with %% correct_bases < FLOAT will be masked. [%f]\n"
            "-e  FLOAT   maximum fraction error bases.  loci with %% error_bases > FLOAT will be masked. [%f]\n\n"
            "Hint: To mask all guide regions, use -c 1.01 -e -0.01 (or anything > 1 and < 0 for -c and -e)\n"
            "      To mask no guide regions, use  -c 0 -e 1e100\n\n"
            "Note: %% is measured as a fraction of guide bases. so\n"
            "      %% correct ranges from [0,1] and\n"
            "      %% error ranges from [0,inf]\n\n",
            cdef, edef);

    fprintf(stderr,
            "jumps.txt has fields:\n\n"
            "<contig> <position> <guide_jump> <correct_jump> <error_jump>\n"
            "See 'align_eval raw' to produce jumps.txt file from SAM alignment\n\n");
    return 1;
}


int main_align_eval_mask(int argc, char ** argv)
{
    float min_frac_correct_bases_def = 0.9;
    float min_frac_correct_bases = min_frac_correct_bases_def;

    float max_frac_error_bases_def = 0.1;
    float max_frac_error_bases = max_frac_error_bases_def;

    char c;
    while ((c = getopt(argc, argv, "c:e:")) >= 0)
    {
        switch(c)
        {
        case 'c': min_frac_correct_bases = atof(optarg); break;
        case 'e': max_frac_error_bases = atof(optarg); break;
        default: return align_eval_mask_usage(min_frac_correct_bases_def,
                                              max_frac_error_bases_def); break;
        }
    }
    
    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return align_eval_mask_usage(min_frac_correct_bases_def,
                                     max_frac_error_bases_def);
    }

    char const* jumps_file = argv[optind];
    char const* mask_file = argv[optind + 1];

    FILE * jumps_fh = open_if_present(jumps_file, "r");
    FILE * mask_fh = open_if_present(mask_file, "w");

    //maintain a boolean value that indicates whether we are in or out
    //of a given mask region.  At each jump bound, we either enter it
    //or exit, based on thresholds
    
    //maintain a running current total of guide, correct, and error.
    char contig[256] = "";
    size_t bound, start_bound = UINT64_MAX;
    int guide_jump, correct_jump, error_jump;
    int guide = 0, correct = 0, error = 0;
    bool in_mask_region = false;
    bool should_be_masked;

    while (! feof(jumps_fh))
    {
        //strcpy(prev_contig, contig);

        fscanf(jumps_fh, "%s\t%zu\t%i\t%i\t%i\n",
               contig, &bound, &guide_jump, &correct_jump, &error_jump);

        guide += guide_jump;
        correct += correct_jump;
        error += error_jump;

        should_be_masked = 
            (guide > 0)
            && ((correct < guide * min_frac_correct_bases)
                || (error > guide * max_frac_error_bases));

        if (in_mask_region)
        {
            //in a mask region.  decide whether to end it.
            if (! should_be_masked)
            {
                //we were previously in a mask region, and this locus
                //should not be masked.  end the mask region
                in_mask_region = false;
                fprintf(mask_fh, "%s\t%Zu\t%Zu\n", contig, start_bound, bound);
            }
        }
        else
        {
            //not in a mask region.  should we start one?
            if (should_be_masked)
            {
                //this locus should be masked.  start it.
                in_mask_region = true;
                start_bound = bound;
            }
        }
    }

    assert(guide == 0);
    assert(correct == 0);
    assert(error == 0);

    close_if_present(mask_fh);
    close_if_present(jumps_fh);

    return 0;
}
