#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

/*
  parse a file of
chrom  locus  value  rest_of_line
chrom  locus  value  rest_of_line
...



produce a file of:

chrom  locus_start locus_end  rest_of_line
chrom  locus_start locus_end  rest_of_line
...

of regions in which the value is >= some threshold given

missing locus values are treated as zero.

the rest of the line can be nothing or something.
if it is something, it must begin with a tab.
i.e., the input and output fields are all tab-separated

 */


int usage()
{
    fprintf(stderr, "Usage\n\n"
            "get_peak_regions threshold_value locus_values.txt peak_regions.txt\n");
    return 0;
}

int main(int argc, char ** argv)
{

    float min_value_for_peak = atof(argv[1]);
    char const* locus_values_file = argv[2];
    char const* peak_regions_file = argv[3];

    FILE * locus_values_fh = fopen(locus_values_file, "r");
    if (locus_values_fh == NULL)
    {
        fprintf(stderr, "Error: Couldn't open input locus values file %s\n",
                locus_values_file);
        exit(1);
    }

    FILE * peak_regions_fh = fopen(peak_regions_file, "w");
    if (peak_regions_fh == NULL)
    {
        fprintf(stderr, "Error: Couldn't open output peak regions file %s\n",
                peak_regions_file);
        exit(1);
    }
    
    bool in_peak = false;
    char cur_chrom[1000];
    char peak_chrom[1000] = "";
    char rest_of_line[1000] = "";

    size_t cur_locus, prev_locus = 0, peak_start = 0;
    float value;

    while (! feof(locus_values_fh))
    {
        //parse locus
        prev_locus = cur_locus;
        fscanf(locus_values_fh, "%s\t%zu\t%f%[^\n]+\n", 
               cur_chrom, &cur_locus, &value, rest_of_line);

        if (rest_of_line[0] != '\0' && rest_of_line[0] != '\t')
        {
            fprintf(stderr, "Error: input file format should be:\n"
                    "chrom<tab>locus<tab>value\n"
                    "or\n"
                    "chrom<tab>locus<tab>value<tab>rest_of_line\n");
            exit(1);
        }
        
        if (in_peak)
        {
            //we were in a peak
            if (value >= min_value_for_peak)
            {
                if (cur_locus == prev_locus + 1 &&
                    strcmp(cur_chrom, peak_chrom) == 0)
                {
                    //the continuation of a peak region.  do nothing
                }
                else
                {
                    //not consecutive.  close and print peak
                    fprintf(peak_regions_fh, "%s\t%Zu\t%Zu%s\n", 
                            peak_chrom, peak_start, prev_locus, rest_of_line);

                    peak_start = cur_locus;
                    strcpy(peak_chrom, cur_chrom);
                    in_peak = true;
                }
            }
            else
            {
                //finish a peak, but do not start a new one
                fprintf(peak_regions_fh, "%s\t%Zu\t%Zu%s\n", 
                        peak_chrom, peak_start, prev_locus, rest_of_line);
                in_peak = false;
            }
        }
        else
        {
            //are in a valley
            if (value >= min_value_for_peak)
            {
                //start a peak
                peak_start = cur_locus;
                strcpy(peak_chrom, cur_chrom);
                in_peak = true;
            }
            else
            {
                //were in a valley, still in it.  do nothing
            }
        }
    }
    if (in_peak)
    {
        //still were in an unclosed peak when we hit the end.
        //print this peak using the current locus, since it wasn't closed yet.
        fprintf(peak_regions_fh, "%s\t%Zu\t%Zu%s\n", 
                peak_chrom, peak_start, cur_locus, rest_of_line);
    }
    fclose(locus_values_fh);
    fclose(peak_regions_fh);
}
