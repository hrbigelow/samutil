#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "samutil.h"

int usage()
{
    fprintf(stderr, "\n"
            "samutil\n"
            "Author: Henry Bigelow (hbigelow@amgen.com)\n\n"
            "Usage:\n\n"
            "samutil sort             sort a SAM, tSAM or rSAM file\n"
            "samutil checksort        check sorted order\n"
            "samutil genome2tx        Condense reads-to-genome alignments to reads-to-transcriptome\n"
            "samutil tx2genome        Expand reads-to-transcriptome alignments to reads-to-genome alignments\n"
            "samutil get_tx_sequence  Get transcript sequence from genomic sequence\n"
            "samutil gen_qcal         Generate qcal file for use with 'samutil score'\n"
            "samutil score            Set mapq, XP, XY, XZ, and primary alignment flags from fragment-sorted SAM file\n"
            "samutil gen_header       Generate a projection-order transcript header\n"      
            "samutil truncate         Truncate a SAM file, replacing QNAME, SEQ, QUAL\n"
            "samutil seqindex         Generate a sequence and index from fastq, to recover QNAME, SEQ, QUAL\n"

            "\n\n"
            "Notes:\n\n"
            "Read names may be in Illumina format, Casava 1.8 foramt, or 'numeric' format\n\n"

            "To ensure that read IDS are truly unique within the set provided,\n"
            "the input must be all of one format\n\n"

            "Illumina format is, e.g.:\n"
            "@HWI-ST630:1:1101:1209:2187#0/1\n"
            "@<flowcell>:<lane>:<tile>:<xpos>:<ypos>[ignored_extra]\n"
            "where all fields except flowcell are positive integers\n\n"

            "Casava 1.8 format:\n"
            "@GAII:1:HWI-ST630:1:1101:1209:2187 1:0:0:ACGGTA\n"
            "@<instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos>[space]\\\n"
            "<read number>:<is filtered>:<control number>:<barcode sequence>\n\n"

            "where all fields except instrument-name and flowcell are positive integers\n\n"

            "numeric format is:\n"
            "@<numeric_id>[ignored_extra]\n\n"
            );
    return 1;
}

int main(int argc, char *argv[])
{
    
    if (argc < 2)
    {
        return usage();
    }
    else if (strcmp(argv[1], "sort") == 0)
    {
        return main_sam_sort(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "checksort") == 0)
    {
        return main_sam_checksort(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "genome2tx") == 0)
    {
        fprintf(stderr, "Sorry, not implemented");
        exit(1);
        //return main_genome_to_transcript(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "tx2genome") == 0)
    {
        return main_transcript_to_genome(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "get_tx_sequence") == 0)
    {
        fprintf(stderr, "Not ported yet");
        exit(1);
        // return main_get_tx_sequence(argc - 1, argv + 1);
    }
    // else if (strcmp(argv[1], "score_dist") == 0)
    // {
    //     return main_score_dist(argc - 1, argv + 1);
    // }
    else if (strcmp(argv[1], "score") == 0)
    {
        return main_score_mapq(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "gen_header") == 0)
    {
        fprintf(stderr, "Not ported yet");
        exit(1);
        //return main_generate_projection_header(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "truncate") == 0)
    {
        return main_sam_truncate(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "seqindex") == 0)
    {
        return main_sam_index_fastq(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "rejoin") == 0)
    {
        return main_sam_rejoin(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
}
