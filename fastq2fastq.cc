//Convert fastq to fastq, with a different qual encoding
#include <algorithm>
#include <cstdio>

#include "fastq_tools.h"

int main(int argc, char ** argv){

    if (argc != 3)
    {
        fprintf(stderr, 
                "Usage: fastq2fastq input.fastq output.fastq\n"
                "\n"
                "Converts to other offset fastq (there are only two unique ones)\n"
                "Assumes input fastq is one of these (from Wikipedia):\n"
                
                "\n"
                "  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................\n"
                "  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................\n"
                "  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................\n"
                "  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................\n"
                "  !\"#$%%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\n"
                "  |                         |    |        |                              |                     |\n"
                " 33                        59   64       73                            104                   126\n"
                "\n"
                " S - Sanger        Phred+33,  raw reads typically (0, 40)\n"
                " X - Solexa        Solexa+64, raw reads typically (-5, 40)\n"
                " I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)\n"
                " J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)\n"
                "\n");
        return 0;

    }

    char * input_file = argv[1];
    char * output_file = argv[2];

    FILE * input_fh = fopen(input_file, "r");
    FILE * output_fh = fopen(output_file, "w");

    if (input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open input file %s\n", input_file);
        exit(1);
    }
    if (output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open output file %s\n", output_file);
        exit(1);
    }

    char id[1024];
    char sequence[1024];
    char quality_string[1024];
    char spacer[1024];
    int nfields_read;
    
    int qual_offset;
    bool fastq_is_valid = fastq_file_offset(input_fh, &qual_offset);

    if (! fastq_is_valid)
    {
        fprintf(stderr, "Error: couldn't determine fastq file quality scale\n");
        exit(1);
    }
    else
    {
        fprintf(stderr, "Input file quality offset is %i:\n", qual_offset);
        rewind(input_fh);
        int input_output_diff = qual_offset == 33 ? (33 - 64) : (64 - 33);
        
        char qseq2fastq[256];
        for (int i=0; i < 256; ++i){
            qseq2fastq[i] = char(std::max(0,i - input_output_diff));
        }
        
        while (! feof(input_fh)){
            nfields_read = 
                fscanf(input_fh, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);
            
            if (nfields_read != 4)
            {
                if (nfields_read > 0)
                {
                    fprintf(stderr, "fastq2fastq: found badly formatted input group with %i lines.\n", nfields_read);
                    exit(1);
                }
                break;
            }
            for (char * q = quality_string; *q != 0; ++q)
            { 
                *q = qseq2fastq[int(*q)];
            }
            
            fprintf(output_fh, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);
        }
        fclose(input_fh);
        fclose(output_fh);
        return 0;
    }
}    
