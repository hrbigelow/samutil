#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdlib>
#include <algorithm>

#include "sam_raw_score_aux.h"
#include "sam_buffer.h"
#include "dep/tools.h"

#include <gsl/gsl_statistics_double.h>

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


typedef std::map<size_t, int> TABLE;
typedef std::map<ScorePair, TABLE> SCORE_CPD;

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
    size_t top_raw_score, sec_raw_score, given_score, num_right, num_wrong;

    double prob_of_correct_alignment;

    SCORE_CPD score_cpd;
    SCORE_CPD::iterator cpd_iter;

    char larger_is_better_char;

    char score_tag[32];
    size_t missing_score_default;
    bool larger_is_better;
    
    num_parsed = fscanf(score_calibration_fh, 
                       "score_tag: %s\n"
                       "missing_score_default: %zu\n"
                       "larger_is_better: %c\n", 
                       score_tag,
                       &missing_score_default,
                       &larger_is_better_char);

    if (num_parsed != 3)
    {
        fprintf(stderr, "Error: score calibration file %s doesn't have score_tag, "
                "missing_score_default or larger_is_better fields.\n"
                "Please produce with 'samutil score_dist'\n",
                score_calibration_file);
        exit(1);
    }

    larger_is_better = larger_is_better_char == 'Y';

    while (! feof(score_calibration_fh))
    {
        num_parsed = fscanf(score_calibration_fh, "%zu\t%zu\t%zu\t%zu\t%zu\n",
                            &top_raw_score, &sec_raw_score, &given_score, &num_right, &num_wrong);
        
        if (num_parsed != 5)
        {
            fprintf(stderr, "Error: score calibration file %s doesn't have five fields per line.\n",
                    score_calibration_file);
            exit(1);
        }

        cpd_iter = score_cpd.find(ScorePair(top_raw_score, sec_raw_score));
        if (cpd_iter == score_cpd.end())
        {
            cpd_iter = 
                score_cpd.insert(std::make_pair(ScorePair(top_raw_score, sec_raw_score), TABLE())).first;
        }
        TABLE & table = (*cpd_iter).second;
        prob_of_correct_alignment = 
            static_cast<double>(num_right) /
            static_cast<double>(num_right + num_wrong);

        table[given_score] = error_prob_to_quality(1.0 - prob_of_correct_alignment);
    }
    fclose (score_calibration_fh);

    //key raw_score sum.  value:  histogram of template lengths
    typedef std::map<size_t, size_t, std::greater<size_t> > HISTO; // general histogram
    typedef std::map<size_t, HISTO, std::greater<size_t> > LENGTH_HISTO_BY_SCORE; // key: score, value: length histo

    LENGTH_HISTO_BY_SCORE length_counts;
    HISTO score_totals;

    size_t total_fragments;

    std::vector<double> template_lengths;

    SAM_ORDER sam_order = SAM_RID_POSITION; //want to catch fragments grouped by physical source fragment
    SamLine::numeric_start_fragment_ids = numeric_read_ids;

    bool paired_reads_are_same_stranded = false;
    bool ones_based_pos = true;
    bool allow_absent_seq_qual = true; // why not?

    int new_mapq;
    size_t template_length;

    bool new_fragment;
    bool ignore_sambuffer_bound;

    bool ignore_duplicate_mapped_pairs = true;

    SamBuffer sam_buffer(sam_order,
                         paired_reads_are_same_stranded,
                         ones_based_pos,
                         ignore_duplicate_mapped_pairs);

    char prev_qname[1024] = "";
    bool seen_a_read = false;
    size_t prev_qid = 0;

    std::vector<size_t> fragment_scores;
    std::vector<size_t>::iterator fit;

    while (! feof(unscored_sam_fh))
    {

        NextLine(unscored_sam_fh, sam_buffer, ones_based_pos, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            ScorePair score_pair = 
                FindTopTwoScores(sam_buffer, 0, SIZE_MAX, score_tag, missing_score_default, 
                                 larger_is_better, & fragment_scores);

            fit = max_element(fragment_scores.begin(), fragment_scores.end());

            size_t max_pair_score = *fit;
            PAIRED_READ_SET::const_iterator mpit = sam_buffer.unique_entry_pairs.begin();
            std::advance(mpit, std::distance(fragment_scores.begin(), fit));

            assert((*mpit).first->isize >= 0);
            template_length = static_cast<size_t>((*mpit).first->isize);
            length_counts[max_pair_score][template_length]++;
            score_totals[max_pair_score]++;
            total_fragments++;
            sam_buffer.purge(NULL, NULL, NULL, ignore_sambuffer_bound);
        }
    }
    rewind(unscored_sam_fh);

    //select the top-scoring N% of fragments as a heuristic proxy for 'correctly aligned' fragments
    size_t partial_sum = 0;
    size_t min_num_used_fragments = 
        static_cast<size_t>(floor(fraction_top_scoring_used * static_cast<float>(total_fragments)));

    HISTO::const_iterator sit_end = score_totals.begin();
    LENGTH_HISTO_BY_SCORE::const_iterator lit_end = length_counts.begin();

    while (partial_sum < min_num_used_fragments)
    {
        partial_sum += (*sit_end++).second;
        ++lit_end;
    }

    //compute mean and sd of elite fragment alignments
    std::vector<double> elite_fragment_lengths;
    LENGTH_HISTO_BY_SCORE::const_iterator lit;
    HISTO::const_iterator hit;
    for (lit = length_counts.begin(); lit != lit_end; ++lit)
    {
        for (hit = (*lit).second.begin(); hit != (*lit).second.end(); ++hit)
        {
            std::fill_n(std::back_inserter(elite_fragment_lengths), (*hit).second, 
                        static_cast<double>((*hit).first));
        }
    }

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

    prev_qid = 0;
    strcpy(prev_qname, "");
    seen_a_read = false;

    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, sam_buffer, ones_based_pos, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            score_pair = FindTopTwoScores(sam_buffer,
                                          min_allowed_fragment_length,
                                          max_allowed_fragment_length,
                                          score_tag,
                                          missing_score_default,
                                          larger_is_better,
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
            
            //4. apply table to each fragment, updating its mapq
            size_t si = 0;
            for (pit = sam_buffer.unique_entry_pairs.begin();
                 pit != sam_buffer.unique_entry_pairs.end(); ++pit)
            {
                SamLine * first = const_cast<SamLine *>((*pit).first);
                SamLine * second = const_cast<SamLine *>((*pit).second);
                table_iter = this_score_table.find(fragment_scores[si++]);
                if (table_iter == this_score_table.end())
                {
                    new_mapq = missing_score_default;
                }
                else
                {
                    new_mapq = (*table_iter).second;
                }
                first->mapq = new_mapq;
                second->mapq = new_mapq;
                
            }
            sam_buffer.purge(scored_sam_fh, NULL, NULL, false);
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
