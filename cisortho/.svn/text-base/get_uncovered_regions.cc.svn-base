#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>


int main(int argc, char ** argv)
{

    int min_window_size = atoi(argv[1]);
    char chrom[1000];
    char last_chrom[1000];
    int64_t gap_start = -1, gap_end = -1, locus = -1;
    int depth = -1;

    int64_t last_locus = -1;
    int64_t last_depth = 1;

    while (! feof(stdin))
    {
  
        scanf("%s\t%"PRId64"\t%i\n", chrom, &locus, &depth);

        if (last_depth == 0
            && (depth > 0
                || strcmp(last_chrom, chrom) != 0))
        {
            //ending a gap
            gap_end = last_locus;
        }
        
        if (gap_start != -1 && gap_end != -1)
        {
            //we have a complete gap
            if (gap_end - gap_start + 1 >= min_window_size)
            {
                //have a large enough gap
                printf("%s\t%"PRId64"\t%"PRId64"\n", last_chrom, gap_start, gap_end);
            }
            gap_start = -1;
            gap_end = -1;
        }

        if (depth == 0 
            && (last_depth > 0 
                || strcmp(last_chrom, chrom) != 0))
        {
            //starting a gap
            gap_start = locus;
        }
        
        strcpy(last_chrom, chrom);
        last_locus = locus;
        last_depth = depth;
    }
  
}
