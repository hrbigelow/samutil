#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdlib>
#include <algorithm>

#include "sam_score_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"

#include <gsl/gsl_statistics_double.h>

int score_mapq_usage(size_t fdef, float qdef, float Qdef, size_t ldef, size_t Ldef, size_t mdef)
{
    fprintf(stderr,
            "\n\nUsage:\n\n"
            "samutil score_mapq [OPTIONS] calibration.qcal unscored.fragsort.sam scored.sam\n\n"
            "Options:\n\n"
            "-f     INT     # top-scoring fragments used for fragment-length distribution [%Zu]\n"
            "-q     FLOAT   min quantile [0,0.5] for -f fragments to calc. min fragment size [%1.2f]\n"
            "-Q     FLOAT   max quantile [0.5,1] for -f fragments to calc. max fragment size [%1.2f]\n"
            "-l     INT     min allowed fragment length for paired alignment [%Zu]\n"
            "-L     INT     max allowed fragment length for paired alignment [%Zu]\n"
            "-n     FLAG    if set, assume numeric read ids (start with an integer) sorting\n"
            "-m     INT     minimum mapq score required to avoid merging top alignments [%Zu]\n"
            "\n\n"
            "calibration.qcal: a histogram over the set of alignment categories\n"
            "(top score, 2nd score, given score)\n"
            "tallying number of alignments that are correct or incorrect.\n\n"

            "unscored_fragsort.sam: alignment file sorted by (read id / pair flag) and\n"
            "having alignment score tags given in option -s.\n\n"

            "scored.sam:    output alignment file with mapq field updated.\n\n"
            "mapq will reflect Phred-scaled probability of correct alignment.\n\n",
            fdef, qdef, Qdef, ldef, Ldef, mdef
            );
    return 1;
}


typedef std::map<size_t, int> TABLE;
typedef std::map<std::pair<size_t, size_t>, TABLE> SCORE_CPD;


int main_score_mapq(int argc, char ** argv)
{
    char c;

    size_t fdef = 100000;
    float qdef = 0.05;
    float Qdef = 0.95;
    size_t ldef = 0;
    size_t Ldef = 1000000;
    size_t mdef = 10;

    size_t num_top_fragments_used = fdef;
    float frag_dist_low_quantile = qdef;
    float frag_dist_high_quantile = Qdef;
    size_t min_allowed_fragment_length = ldef;
    size_t max_allowed_fragment_length = Ldef;
    size_t min_mapq_for_not_merging = mdef;

    bool numeric_read_ids = false;

    while ((c = getopt(argc, argv, "f:q:Q:l:L:nm:")) >= 0)
    {
        switch(c)
        {
        case 'f': num_top_fragments_used = static_cast<size_t>(atoi(optarg)); break;
        case 'q': frag_dist_low_quantile = atof(optarg); break;
        case 'Q': frag_dist_high_quantile = atof(optarg); break;
        case 'l': min_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'L': max_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'n': numeric_read_ids = true; break;
        case 'm': min_mapq_for_not_merging = static_cast<size_t>(atoi(optarg)); break;
        default: return score_mapq_usage(fdef, qdef, Qdef, ldef, Ldef, mdef); break;
        }
    }

    if (argc != optind + 3)
    {
        return score_mapq_usage(fdef, qdef, Qdef, ldef, Ldef, mdef);
    }

    char * score_calibration_file = argv[optind];
    char * unscored_sam_file = argv[optind + 1];
    char * scored_sam_file = argv[optind + 2];

    FILE * score_calibration_fh = open_or_die(score_calibration_file, "r", "Input score calibration file");
    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * scored_sam_fh = open_or_die(scored_sam_file, "w", "Output scored sam file");

    //parse score calibration
    int num_parsed;
    size_t top_raw_score, sec_raw_score, given_raw_score, num_right, num_wrong;

    double prob_of_correct_alignment;

    SCORE_CPD score_cpd;
    SCORE_CPD::iterator cpd_iter;

    char larger_score_better_char;

    char score_tag[32];
    size_t missing_default_score;


    num_parsed = fscanf(score_calibration_fh, 
                       "score_tag: %s\n"
                       "missing_default_score: %zu\n"
                       "larger_score_better: %c\n"
                        score_tag,
                        &missing_default_score,
                        &larger_score_better_char);

    if (num_parsed != 4)
    {
        fprintf(stderr, "Error: score calibration file %s doesn't have score_tag, "
                "missing_default_score, or larger_is_better fields.\n"
                "Please produce with 'samutil score_dist'\n",
                score_calibration_file);
        exit(1);
    }

    bool larger_score_better = larger_score_better_char == 'Y';

    while (! feof(score_calibration_fh))
    {
        num_parsed = fscanf(score_calibration_fh, "%zu\t%zu\t%zu\t%zu\t%zu\n",
                            &top_raw_score,
                            &sec_raw_score,
                            &given_raw_score,
                            &num_right, &num_wrong);
        
        if (num_parsed != 5)
        {
            fprintf(stderr, "Error: score calibration file %s doesn't have five fields per line.\n",
                    score_calibration_file);
            exit(1);
        }

        cpd_iter = score_cpd.find(std::make_pair(top_raw_score, sec_raw_score));

        if (cpd_iter == score_cpd.end())
        {
            cpd_iter = 
                score_cpd.insert(std::make_pair(ScorePair(top_raw_score, sec_raw_score), TABLE())).first;
        }
        TABLE & table = (*cpd_iter).second;
        prob_of_correct_alignment = 
            static_cast<double>(num_right) /
            static_cast<double>(num_right + num_wrong);

        table[given_raw_score] = error_prob_to_quality(1.0 - prob_of_correct_alignment);
    }
    fclose (score_calibration_fh);

    //use as the 'less' comparator to consider better scores to be 'less'
    //Places the best scores at the beginning.
    RawScoreBetter score_order(larger_score_better);

    //key raw_score sum.  value:  histogram of template lengths

    std::vector<double> template_lengths;

    SamLine::SetGlobalFlags(numeric_read_ids);

    SamOrder sam_order(SAM_RID_POSITION, "NONE");
    sam_order.AddHeaderContigStats(unscored_sam_fh);

    bool paired_reads_are_same_stranded = false;
    bool allow_absent_seq_qual = true; // why not?

    int new_mapq;

    bool new_fragment;
    bool ignore_sambuffer_bound;

    SamBuffer tally_buffer(&sam_order, paired_reads_are_same_stranded);

    //tally fragment length and raw score statistics
    size_t min_fragment_length;
    size_t max_fragment_length;

    QuantileFragmentEstimate(min_allowed_fragment_length,
                             max_allowed_fragment_length,
                             frag_dist_low_quantile,
                             frag_dist_high_quantile,
                             &unscored_sam_fh, 
                             score_order, 
                             tally_buffer, 
                             score_tag,
                             missing_default_score,
                             num_top_fragments_used,
                             &min_fragment_length,
                             &max_fragment_length);

    fprintf(stderr, 
            "Fragment length number top-scoring alignments used: %Zu\n"
            "Fragment length quantiles used: %f - %f\n"
            "Fragment length hard limits: %Zu - %Zu\n"
            "Fragment length range determined as: %Zu - %Zu\n",
            num_top_fragments_used,
            frag_dist_low_quantile, frag_dist_high_quantile,
            min_allowed_fragment_length, max_allowed_fragment_length,
            min_fragment_length, max_fragment_length);

    fflush(stderr);

    PAIRED_READ_SET::iterator pit;
    TABLE::const_iterator table_iter;
    TABLE const* table_ptr;
    TABLE empty_table;

    PrintSAMHeader(&unscored_sam_fh, scored_sam_fh);
    fflush(scored_sam_fh);

    size_t prev_qid = 0;
    char prev_qname[1024] = "";
    bool seen_a_read = false;

    SamBuffer cal_buffer(&sam_order, paired_reads_are_same_stranded);
    
    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, cal_buffer, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            SetScoreFields(cal_buffer,
                           min_fragment_length,
                           max_fragment_length,
                           raw_score_tag,
                           missing_default_score,
                           larger_raw_score_is_better,
                           score_cpd);
            
            cal_buffer.purge(scored_sam_fh, NULL, NULL, false);
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
