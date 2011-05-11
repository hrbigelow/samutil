#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdlib>
#include <algorithm>

#include "sam_score_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_uint.h>

int score_mapq_usage(float fdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil score_mapq [OPTIONS] raw_score_calibration.txt unscored.fragsort.sam scored.sam\n\n"
            "Options:\n\n"
            "-f      FLOAT    fraction of top-scoring alignments to use to determine fragment length distribution [%f]\n"
            "-n      FLAG     if set, assume numeric read ids (start with an integer) sorting [false]\n"
            "\n\n"
            "raw_score_calibration.txt: a probability table\n"
            "P(correct alignment, raw score | top raw score, 2nd raw score)\n"
            "Format:  top_raw_score 2nd_raw_score given_score probability_of_correct_alignment\n"
            "unscored.sam:  alignment file sorted by fragment identity (read id / pair flag)\n"
            "having raw score tag 'AS:i'.\n"
            "scored.sam:    output alignment file with mapq field updated.\n\n"
            "mapq will reflect Phred-scaled probability of correct alignment.\n"
            "In case given_score is below the top two scores, mapq is assigned zero\n",
            fdef
            );
    return 1;
}


typedef std::map<SpaceScore, int> TABLE;
typedef std::map<ScorePair, TABLE> SCORE_CPD;


struct ScoreOrder
{
    bool increasing;
    ScoreOrder(bool _inc = true) : increasing(_inc) { }
    bool operator()(size_t a, size_t b)
    {
        return this->increasing ? a < b : b < a;
    }
};


int main_score_mapq(int argc, char ** argv)
{
    char c;

    float fdef = 0.05;
    float fraction_top_scoring_used = fdef;

    bool numeric_read_ids = false;

    while ((c = getopt(argc, argv, "f:n")) >= 0)
    {
        switch(c)
        {
        case 'f': fraction_top_scoring_used = atof(optarg); break;
        case 'n': numeric_read_ids = true; break;
        default: return score_mapq_usage(fdef); break;
        }
    }

    if (argc != optind + 3)
    {
        return score_mapq_usage(fdef);
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
    char larger_space_better_char;

    char score_tag[32];
    size_t missing_default_score;


    num_parsed = fscanf(score_calibration_fh, 
                       "score_tag: %s\n"
                       "missing_default_score: %zu\n"
                       "larger_score_better: %c\n"
                       "larger_space_better: %c\n", 
                        score_tag,
                        &missing_default_score,
                        &larger_score_better_char,
                        &larger_space_better_char);

    if (num_parsed != 4)
    {
        fprintf(stderr, "Error: score calibration file %s doesn't have score_tag, "
                "missing_default_score, larger_is_better, or larger_space_better fields.\n"
                "Please produce with 'samutil score_dist'\n",
                score_calibration_file);
        exit(1);
    }

    bool larger_score_better = larger_score_better_char == 'Y';
    bool larger_space_better = larger_space_better_char == 'Y';

    char top_raw_space;
    char sec_raw_space;
    char given_raw_space;

    while (! feof(score_calibration_fh))
    {
        num_parsed = fscanf(score_calibration_fh, "%zu\t%c\t%zu\t%c\t%zu\t%c\t%zu\t%zu\n",
                            &top_raw_score, &top_raw_space, 
                            &sec_raw_score, &sec_raw_space,
                            &given_raw_score, &given_raw_space,
                            &num_right, &num_wrong);
        
        if (num_parsed != 8)
        {
            fprintf(stderr, "Error: score calibration file %s doesn't have eight fields per line.\n",
                    score_calibration_file);
            exit(1);
        }

        AlignmentSpace top_space = ParseAlignmentSpace(top_raw_space);
        AlignmentSpace sec_space = ParseAlignmentSpace(sec_raw_space);
        AlignmentSpace given_space = ParseAlignmentSpace(given_raw_space);

        SpaceScore top_spacescore = SpaceScore(top_raw_score, top_space);
        SpaceScore sec_spacescore = SpaceScore(sec_raw_score, sec_space);
        SpaceScore given_spacescore = SpaceScore(given_raw_score, given_space);

        cpd_iter = score_cpd.find(ScorePair(top_spacescore, sec_spacescore));

        if (cpd_iter == score_cpd.end())
        {
            cpd_iter = 
                score_cpd.insert(std::make_pair(ScorePair(top_spacescore, sec_spacescore), TABLE())).first;
        }
        TABLE & table = (*cpd_iter).second;
        prob_of_correct_alignment = 
            static_cast<double>(num_right) /
            static_cast<double>(num_right + num_wrong);

        table[given_spacescore] = error_prob_to_quality(1.0 - prob_of_correct_alignment);
    }
    fclose (score_calibration_fh);

    //use as the 'less' comparator to consider better scores to be 'less'
    //Places the best scores at the beginning.
    SpaceScore_better score_best_is_less(larger_score_better, larger_space_better);

    //key raw_score sum.  value:  histogram of template lengths

    std::vector<double> template_lengths;

    SAM_ORDER sam_order = SAM_RID_POSITION; //want to catch fragments grouped by physical source fragment
    SamLine::numeric_start_fragment_ids = numeric_read_ids;

    bool paired_reads_are_same_stranded = false;
    bool allow_absent_seq_qual = true; // why not?

    int new_mapq;

    bool new_fragment;
    bool ignore_sambuffer_bound;

    bool ignore_duplicate_mapped_pairs = true;

    

    ScoreVec fragment_scores;
    ScoreVec::iterator fit;

    SamBuffer tally_buffer(sam_order,
                           paired_reads_are_same_stranded,
                           ignore_duplicate_mapped_pairs);

    //tally fragment length and raw score statistics
    std::vector<double> elite_fragment_lengths =
        TallyFragmentLengths(&unscored_sam_fh, 
                             score_best_is_less, 
                             tally_buffer, 
                             score_tag,
                             missing_default_score,
                             fraction_top_scoring_used);

    double fragment_mean = gsl_stats_mean(elite_fragment_lengths.data(), 1,
                                          elite_fragment_lengths.size());

    double fragment_sd = gsl_stats_sd_m(elite_fragment_lengths.data(), 1,
                                        elite_fragment_lengths.size(), fragment_mean);

    size_t min_allowed_fragment_length = static_cast<size_t>(floor(fragment_mean - (2.0 * fragment_sd)));
    size_t max_allowed_fragment_length = static_cast<size_t>(floor(fragment_mean + (2.0 * fragment_sd)));

    PAIRED_READ_SET::iterator pit;
    ScorePair score_pair;
    TABLE::const_iterator table_iter;
    TABLE const* table_ptr;
    TABLE empty_table;

    PrintSAMHeader(&unscored_sam_fh, scored_sam_fh);
    fflush(scored_sam_fh);


    size_t prev_qid = 0;
    char prev_qname[1024] = "";
    bool seen_a_read = false;


    SamBuffer cal_buffer(sam_order,
                         paired_reads_are_same_stranded,
                         ignore_duplicate_mapped_pairs);
    

    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, cal_buffer, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            score_pair = FindTopTwoScores(cal_buffer,
                                          min_allowed_fragment_length,
                                          max_allowed_fragment_length,
                                          score_tag,
                                          missing_default_score,
                                          score_best_is_less,
                                          &fragment_scores);
            
            cpd_iter = score_cpd.find(score_pair);
            
            if (cpd_iter == score_cpd.end())
            {
                //no information about how to calculate mapq from raw scores.
                //resort to an empty table, which will trigger using the low default
                table_ptr = &empty_table;
            }
            else
            {
                table_ptr = &(*cpd_iter).second;
            }
            TABLE const& this_score_table = *table_ptr;
            
            //4. apply table to each fragment, updating its mapq and HI tag field
            size_t si = 0;

            //use std::greater because mapq is defined as 'higher-is-better', and
            //rank is defined as 'lower-is-better'
            std::multimap<int, size_t, std::greater<int> > mapqs;
            std::multimap<int, size_t, std::greater<int> >::const_iterator mit;

            size_t N = cal_buffer.unique_entry_pairs.size();

            uint * ranks = new uint[N];
            gsl_permutation * perm = gsl_permutation_calloc(N);

            for (pit = cal_buffer.unique_entry_pairs.begin();
                 pit != cal_buffer.unique_entry_pairs.end(); ++pit)
            {
                SamLine * first = const_cast<SamLine *>((*pit).first);
                SamLine * second = const_cast<SamLine *>((*pit).second);
                table_iter = this_score_table.find(fragment_scores[si]);
                if (table_iter == this_score_table.end())
                {
                    // we're in a state of complete ignorance.  This shouldn't happen.
                    new_mapq = 0;
                    assert(false);
                }
                else
                {
                    new_mapq = (*table_iter).second;
                }
                first->mapq = new_mapq;
                second->mapq = new_mapq;
                
                mapqs.insert(std::make_pair(first->mapq, si));
                ranks[si] = static_cast<uint>(si);
                ++si;
            }

            for (mit = mapqs.begin(), si = 0; mit != mapqs.end(); ++mit)
            {
                perm->data[(*mit).second] = si;
                ++si;
            }

            gsl_permute_uint(perm->data, ranks, 1, N);

            char tag_template[2000];
            si = 0;
            for (pit = cal_buffer.unique_entry_pairs.begin();
                 pit != cal_buffer.unique_entry_pairs.end(); ++pit)
            {
                sprintf(tag_template, "HI:i:%u", ranks[si] + 1);
                SamLine * first = const_cast<SamLine *>((*pit).first);
                SamLine * second = const_cast<SamLine *>((*pit).second);
                first->add_tag(tag_template);
                second->add_tag(tag_template);
                if (ranks[si] == 0 && first->mapped_in_proper_pair())
                {
                    //alignment is primary
                    first->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
                    second->flag &= ~SamFlags::ALIGNMENT_NOT_PRIMARY;
                }
                else
                {
                    first->flag |= SamFlags::ALIGNMENT_NOT_PRIMARY;
                    second->flag |= SamFlags::ALIGNMENT_NOT_PRIMARY;
                }
                ++si;
            }   
         
            cal_buffer.purge(scored_sam_fh, NULL, NULL, false);
            delete ranks;
            gsl_permutation_free(perm);
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
