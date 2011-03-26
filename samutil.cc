#include <cstdio>
#include <cstring>

#include "samutil.h"

int usage()
{
    fprintf(stderr, "\n"
            "samutil\n"
            "Author: Henry Bigelow (hbigelow@amgen.com)\n\n"
            "Usage:\n\n"
            "samutil genome2tx        Condense reads-to-genome alignments to reads-to-transcriptome\n"
            "samutil tx2genome        Expand reads-to-transcriptome alignments to reads-to-genome alignments\n"
            "samutil get_tx_sequence  Get transcript sequence from genomic sequence\n"
            "samutil score_dist       Generate raw score distribution statistics for 'score_mapq'\n"
            "samutil score_mapq       Score alignment mapping quality from fragment-sorted SAM file\n"
            "\n"
            );
    return 1;
}

int main(int argc, char *argv[])
{
    
    if (argc < 2)
    {
        return usage();
    }
    else if (strcmp(argv[1], "genome2tx") == 0)
    {
        return main_genome_to_transcript(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "tx2genome") == 0)
    {
        return main_transcript_to_genome(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "get_tx_sequence") == 0)
    {
        return main_get_tx_sequence(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "score_dist") == 0)
    {
        return main_score_dist(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "score_mapq") == 0)
    {
        return main_score_mapq(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
}
