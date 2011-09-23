#include "align_eval_raw.h"
#include "align_eval_mask.h"
#include "align_eval_coverage.h"
#include "align_eval_stats.h"

int main_usage()
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "align_eval raw [OPTIONS] alignment_sorted.sam jumps.txt align_stats.cumul.txt align_stats.{oplen,fragsize}.{full,by_half}.txt \n"
            "align_eval mask [OPTIONS] jumps.txt mask.txt\n"
            "align_eval coverage [OPTIONS] jumps.txt coverage.txt\n"
            "align_eval stats [OPTIONS] jumps.txt stats.txt dist.txt\n\n"
            );

    return 1;
}


int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        return main_usage();
    }
    else if (strcmp(argv[1], "raw") == 0)
    {
        return main_align_eval_raw(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "mask") == 0)
    {
        return main_align_eval_mask(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "coverage") == 0)
    {
        return main_align_eval_coverage(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "stats") == 0)
    {
        return main_align_eval_stats(argc - 1, argv + 1);
    }
    else
    {
        return main_usage();
    }
}
