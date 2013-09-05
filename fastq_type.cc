#include "fastq_tools.h"


int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        fprintf(stderr,
                "Usage: fastq_type input.fastq\n\n"
                "Scans input fastq, outputs:\n"
                "<valid|invalid>\t<quality_code_offset>\t<file_types>\n"
                "e.g.:\n"
                "valid\t64\t-XI-\n"
                "<file_types> is of format [S-][X-][I-][J-]\n"
                "A letter indicates that all QUALs fit within the valid range.\n"
                "A dash indicates that one or more QUAL characters fell outside the range (as defined below).\n"
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
                " J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)\n\n"
                );
        return 0;
    }

    FILE * fh = fopen(argv[1], "r");

    char minc, maxc;
    fastq_extreme_chars(fh, & minc, & maxc);
    
    fprintf(stdout, "%c-%c\n", minc, maxc);

    fclose(fh);

    return 0;
}
