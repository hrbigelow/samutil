#include <map>
#include <cstdio>
#include <cstdlib>

#include "sam_raw_score_aux.h"
#include "dep/tools.h"
#include "sam_helper.h"

/*
generating a JPD for raw scores
1. simulate and align paired reads with favorite aligner
2. sort by read_id, read_pair
3. parse all alignments of each given fragment, one by one.
   for each:
   a.  determine the top two raw fragment scores (sum of raw read alignment scores)
   b.  using a set base accuracy threshold, call each alignment 'correct' or 'incorrect',
   and update a histogram of (top_raw_score, 2nd_raw_score, given_score, correct_count, total_count)
   (by convention, if there is only one alignment, use 0 for the 2nd_raw_score)

4. output this table.
 */


int raw_score_dist_usage(float fdef, size_t ndef, size_t xdef, char const* sdef, size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil score_dist [OPTIONS] unscored_fragsort.sam raw_score_calibration.txt\n\n"
            "Options:\n\n"
            "-f      FLOAT    fraction of correctly aligned bases to call an alignment 'correct' [%f]\n"
            "-n      INT      minimum allowed fragment length for paired alignment [%Zu]\n"
            "-x      INT      maximum allowed fragment length for paired alignment [%Zu]\n"
            "-s      STRING   SAM tag representing the alignment score.  Must be an integer valued tag [%s]\n"
            "-l      FLAG     if present, consider that a lower alignment score is better [false]\n"
            "-m      INT      default alignment score assumed for SAM entries with missing tag [%Zu]\n"
            "\n\n"
            "raw_score_calibration.txt: a 4D histogram with\n"
            "(top raw score, 2nd raw score, given score, is_correct)\n"
            "unscored.fragsort.sam:  alignment file sorted by fragment identity (read id / pair flag)\n"
            "having raw score tag 'AS:i'.\n",
            fdef, ndef, xdef, sdef, mdef
            );
    return 1;
}


struct ScoreCategory
{
    size_t top_score;
    size_t sec_score;
    size_t given_score;
    ScoreCategory(size_t _t, size_t _s, size_t _g) : 
        top_score(_t), sec_score(_s), given_score(_g) { }
    ScoreCategory() : 
        top_score(0), sec_score(0), given_score(0) { }
    bool operator<(ScoreCategory const& b) const
    {
        return this->top_score < b.top_score
            || (this->top_score == b.top_score
                && this->sec_score < b.sec_score
                || (this->sec_score == b.sec_score
                    && this->given_score < b.given_score
                    ));
    }
};


typedef std::map<ScoreCategory, std::pair<size_t, size_t> > RAW_HISTO;

int main_score_dist(int argc, char ** argv)
{
    char c;
    float fdef = 1.0;
    size_t ndef = 100;
    size_t xdef = 300;
    char const* sdef = "AS";
    size_t mdef = 0;

    float min_correctly_aligned_bases = fdef;
    size_t min_allowed_fragment_length = ndef;
    size_t max_allowed_fragment_length = xdef;
    char const* score_tag = sdef;
    size_t missing_score_default = mdef;
    bool larger_is_better = true;
    

    while ((c = getopt(argc, argv, "f:n:x:s:lm:")) >= 0)
    {
        switch(c)
        {
        case 'f': min_correctly_aligned_bases = atof(optarg); break;
        case 'n': min_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'x': max_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 's': score_tag = optarg; break;
        case 'l': larger_is_better = false; break;
        case 'm': missing_score_default = static_cast<size_t>(atoi(optarg)); break; 
        default: return raw_score_dist_usage(fdef, ndef, xdef, sdef, mdef); break;
        }
    }

    if (argc != optind + 2)
    {
        return raw_score_dist_usage(fdef, ndef, xdef, sdef, mdef);
    }

    char * unscored_sam_file = argv[optind];
    char * score_calibration_file = argv[optind + 1];

    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * score_calibration_fh = open_or_die(score_calibration_file, "w", "Output score calibration file");

    SAM_ORDER sam_order = SAM_RID_POSITION; //want to catch fragments
    SamLine::numeric_start_fragment_ids = true; // don't assume we have guide data.

    //alignments from different physical fragments
    bool paired_reads_are_same_stranded = false;
    bool ones_based_pos = true;
    bool allow_absent_seq_qual = true; // why not?
    bool ignore_duplicate_mapped_pairs = true; //anticipate a mixture of alignments from different
    //sources, and take the unique subset of them

    int missing_second_raw_score = 0;
    int new_mapq;
    bool new_fragment;

    bool ignore_sambuffer_bound;

    SamBuffer sam_buffer(sam_order,
                         paired_reads_are_same_stranded,
                         ones_based_pos,
                         ignore_duplicate_mapped_pairs);

    char prev_qname[1024] = "";

    std::vector<size_t> fragment_scores;

    ScorePair score_pair;

    RAW_HISTO score_cat_histo;

    read_coords first_guide_coords;
    read_coords second_guide_coords;

    SamLine *samline;
    bool seen_a_read = false;
    
    size_t prev_qid = 0;

    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, sam_buffer, ones_based_pos, allow_absent_seq_qual,
                 &new_fragment, &ignore_sambuffer_bound,
                 &seen_a_read, prev_qname, &prev_qid);

        if (new_fragment)
        {
            //find scores of top two fragments.
            score_pair = FindTopTwoScores(sam_buffer,
                                          min_allowed_fragment_length,
                                          max_allowed_fragment_length,
                                          score_tag, missing_score_default,
                                          larger_is_better,
                                          &fragment_scores);

            //then, for each fragment, evaluate whether it is correctly aligned,
            //and update the histogram

            PAIRED_READ_SET::const_iterator pit = sam_buffer.unique_entry_pairs.begin();

            CigarFromSimSAMLine((*pit).first->qname, 
                                (*pit).first->first_read_in_pair(),
                                ones_based_pos,
                                &first_guide_coords);

            CigarFromSimSAMLine((*pit).second->qname, 
                                (*pit).second->first_read_in_pair(),
                                ones_based_pos,
                                &second_guide_coords);
            
            Cigar::CIGAR_VEC first_guide_cigar = 
                Cigar::FromString(first_guide_coords.cigar, first_guide_coords.position);

            std::multiset<std::pair<size_t, size_t> > first_guide_cigar_index =
                Cigar::ComputeOffsets(first_guide_cigar);

            Cigar::CIGAR_VEC second_guide_cigar = 
                Cigar::FromString(second_guide_coords.cigar, second_guide_coords.position);

            std::multiset<std::pair<size_t, size_t> > second_guide_cigar_index =
                Cigar::ComputeOffsets(second_guide_cigar);

            size_t frag_index = 0;
            for (pit = sam_buffer.unique_entry_pairs.begin();
                 pit != sam_buffer.unique_entry_pairs.end(); ++pit)
            {
                SamLine const* first = (*pit).first;
                SamLine const* second = (*pit).second;

                size_t first_num_correct_bases = 
                    CountCorrectBases(first, first_guide_coords,
                                      first_guide_cigar,
                                      first_guide_cigar_index);
                
                size_t second_num_correct_bases = 
                    CountCorrectBases(second, second_guide_coords,
                                      second_guide_cigar,
                                      second_guide_cigar_index);
                
                size_t num_correctly_aligned_bases = 
                    first_num_correct_bases
                    + second_num_correct_bases;

                bool is_correct = num_correctly_aligned_bases >= min_correctly_aligned_bases;

                ScoreCategory score_cat(score_pair.top_score, score_pair.sec_score,
                                        fragment_scores[frag_index++]);

                if (score_cat_histo.find(score_cat) == score_cat_histo.end())
                {
                    score_cat_histo[score_cat] = std::make_pair(0, 0);
                }

                if (is_correct)
                {
                    score_cat_histo[score_cat].first++;
                }
                else
                {
                    score_cat_histo[score_cat].second++;
                }
            }

            sam_buffer.purge(NULL, NULL, NULL, ignore_sambuffer_bound);
        }
    }

    RAW_HISTO::const_iterator rit;
    fprintf(score_calibration_fh, "score_tag: %s\n", score_tag);
    fprintf(score_calibration_fh, "missing_score_default: %zu\n", missing_score_default);
    fprintf(score_calibration_fh, "larger_is_better: %c\n", (larger_is_better ? 'Y' : 'N'));

    for (rit = score_cat_histo.begin(); rit != score_cat_histo.end(); ++rit)
    {
        ScoreCategory const& cat = (*rit).first;
        size_t num_right = (*rit).second.first;
        size_t num_wrong = (*rit).second.second;

        fprintf(score_calibration_fh, "%Zu\t%Zu\t%Zu\t%Zu\t%Zu\n",
                cat.top_score, cat.sec_score, cat.given_score, num_right, num_wrong);
    }

    fclose(score_calibration_fh);
    fclose(unscored_sam_fh);

    return 0;
}
