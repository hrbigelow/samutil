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
            "samutil score_mapq [OPTIONS] calibration.qcal unscored.rsort.sam scored.rsort.sam\n\n"
            "Options:\n\n"
            "-l     INT     min allowed fragment length for paired alignment [%Zu]\n"
            "-L     INT     max allowed fragment length for paired alignment [%Zu]\n"
            "-e     FLAG    if present, use equivalency mapq. if absent, use traditional [absent]\n"
            "\n\n"
            "calibration.qcal: a histogram over the set of alignment categories\n"
            "(top score, 2nd score, given score)\n"
            "tallying number of alignments that are correct or incorrect.\n\n"

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
    // bool equivalency_mapq = false;

    while ((c = getopt(argc, argv, "l:L:")) >= 0)
    {
        switch(c)
        {
        case 'l': min_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'L': max_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        // case 'e': equivalency_mapq = true; break;
        default: return score_mapq_usage(ldef, Ldef); break;
        }
    }

    if (argc != optind + 3)
    {
        return score_mapq_usage(ldef, Ldef);
    }

    char * score_calibration_file = argv[optind];
    char * unscored_sam_file = argv[optind + 1];
    char * scored_sam_file = argv[optind + 2];

    FILE * score_calibration_fh = open_or_die(score_calibration_file, "r", "Input score calibration file");
    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * scored_sam_fh = open_or_die(scored_sam_file, "w", "Output scored sam file");

    char raw_score_tag[32];
    bool larger_score_better;

    int * score_cpd;

    CollapseScoreTriplet score_triplet =
        ParseScoreCalibration(score_calibration_fh,
                              raw_score_tag,
                              &larger_score_better,
                              score_cpd);

    RAW_SCORE_T max_valid_score = score_triplet.max_score - 1;

    fclose(score_calibration_fh);

    //use as the 'less' comparator to consider better scores to be 'less'
    //Places the best scores at the beginning.
    RawScoreBetter score_order(larger_score_better);

    //key raw_score sum.  value:  histogram of template lengths

    std::vector<double> template_lengths;

    SamOrder sam_order(SAM_RID_POSITION, "NONE");
    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(unscored_sam_fh);
    sam_order.AddHeaderContigStats(unscored_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt);

    bool paired_reads_are_same_stranded = false;
    bool allow_absent_seq_qual = true; // why not?

    bool new_fragment;

    SamBuffer tally_buffer(&sam_order, paired_reads_are_same_stranded);

    PAIRED_READ_SET::iterator pit;

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
            SetScoreFields(cal_buffer,
                           min_fragment_length,
                           max_fragment_length,
                           raw_score_tag,
                           max_valid_score,
                           larger_score_better,
                           score_cpd,
                           score_triplet);
            
            cal_buffer.purge(scored_sam_fh, NULL, NULL, low_bound);
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
