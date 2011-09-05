#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdlib>
#include <algorithm>

#include "sam_score_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"

#include <gsl/gsl_statistics_double.h>

int score_mapq_usage(size_t ldef, size_t Ldef)
{
    fprintf(stderr,
            "\n\nUsage:\n\n"
            "samutil score [OPTIONS] calibration.qcal contig_space.txt unscored.rsort.sam scored.rsort.sam\n\n"
            "Options:\n\n"
            "-l     INT     min allowed fragment length for paired alignment [%Zu]\n"
            "-L     INT     max allowed fragment length for paired alignment [%Zu]\n"
            "-r     FLAG    If present, print in rSAM format. Otherwise, print traditional SAM [false]\n"
            "\n\n"
            "calibration.qcal: a histogram over the set of alignment categories\n"
            "(top score, 2nd score, given score)\n"
            "tallying number of alignments that are correct or incorrect.\n\n"

            "contig_space.txt: defines contig alignment space for each contig in 'unscored.rsort.sam' file\n"
            "Format:\n\n"

            "# comments ...\n"
            "TG\n"
            "chr1<tab>G\n"
            "chr2<tab>G\n"
            "...\n"
            "ENST000005299<tab>T\n"
            "ENST899998402<tab>T\n"
            "...\n\n"

            "unscored_rsort.sam: alignment file sorted by (read id / pair flag) and\n"
            "having alignment score tags given in option -s.\n\n"

            "scored.rsort.sam:    output alignment file sorted by (read id / pair flag)\n"
            "with mapq field updated.\n"
            "mapq will reflect Phred-scaled probability of correct alignment.\n\n",
            ldef, Ldef
            );
    return 1;
}



int main_score_mapq(int argc, char ** argv)
{
    char c;

    size_t ldef = 0;
    size_t Ldef = 1000;

    size_t min_fragment_length = ldef;
    size_t max_fragment_length = Ldef;
    bool print_rsam = false;

    // bool equivalency_mapq = false;

    while ((c = getopt(argc, argv, "l:L:")) >= 0)
    {
        switch(c)
        {
        case 'l': min_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'L': max_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'r': print_rsam = true; break;
        // case 'e': equivalency_mapq = true; break;
        default: return score_mapq_usage(ldef, Ldef); break;
        }
    }

    if (argc != optind + 4)
    {
        return score_mapq_usage(ldef, Ldef);
    }

    char * score_calibration_file = argv[optind];
    char * contig_space_file = argv[optind + 1];
    char * unscored_sam_file = argv[optind + 2];
    char * scored_sam_file = argv[optind + 3];

    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * scored_sam_fh = open_or_die(scored_sam_file, "w", "Output scored sam file");

    FragmentScore fragment_scoring(min_fragment_length, max_fragment_length);
    fragment_scoring.init(score_calibration_file, contig_space_file);


    SamOrder sam_order(SAM_RID_POSITION, "NONE");
    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(unscored_sam_fh);
    sam_order.AddHeaderContigStats(unscored_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt);

    bool paired_reads_are_same_stranded = false;
    bool allow_absent_seq_qual = true; // why not?

    bool new_fragment;

    SamBuffer tally_buffer(&sam_order, paired_reads_are_same_stranded);

    PrintSAMHeader(&unscored_sam_fh, scored_sam_fh);
    fflush(scored_sam_fh);

    size_t prev_qid = 0;
    char prev_qname[1024] = "";
    bool seen_a_read = false;

    SamBuffer cal_buffer(&sam_order, paired_reads_are_same_stranded);

    SamLine * low_bound = NULL;

    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, cal_buffer, allow_absent_seq_qual,
                 &new_fragment, &seen_a_read, prev_qname, &prev_qid, &low_bound);

        if (new_fragment)
        {
            set_score_fields(cal_buffer, fragment_scoring);
            cal_buffer.purge(scored_sam_fh, NULL, NULL, print_rsam, low_bound);
        }
        else
        {
            //not a new fragment.  just continue loading.
        }
    }

    fclose(unscored_sam_fh);
    fclose(scored_sam_fh);
    return 0;
}
