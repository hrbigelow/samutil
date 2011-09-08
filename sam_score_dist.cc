#include <map>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "sam_score_aux.h"
#include "dep/tools.h"
#include "sam_helper.h"
#include "sam_order.h"

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


int score_dist_usage(size_t fdef, float qdef, float Qdef, size_t ldef, size_t Ldef, 
                     float cdef, char const* sdef, size_t mdef)
{
    fprintf(stderr,
            "\n\nUsage:\n\n"
            "samutil score_dist [OPTIONS] unscored.fragsort.sam calibration.qcal\n\n"
            "Options:\n\n"
            "-f     INT     # top-scoring fragments used for fragment-length distribution [%Zu]\n"
            "-q     FLOAT   min quantile [0,0.5] for -f fragments to calc. min fragment size [%1.2f]\n"
            "-Q     FLOAT   max quantile [0.5,1] for -f fragments to calc. max fragment size [%1.2f]\n"
            "-l     INT     min allowed fragment length for paired alignment [%Zu]\n"
            "-L     INT     max allowed fragment length for paired alignment [%Zu]\n"
            "-c     FLOAT   fraction of correctly aligned bases to call an alignment 'correct' [%1.2f]\n"
            "-s     STRING  SAM tag representing the alignment score.  Must be \":i\" tag [%s]\n"
            "-i     FLAG    if present, consider a larger alignment score better [false]\n"
            "-m     INT     default alignment score assumed for SAM entries with missing tag [%Zu]\n"
            "-e     FLAG    if set, alignment is 'correct' if given raw score is same as that of correct alignment [false]\n"
            "-y     STRING   expected read layout. If parsing traditional SAM, this is required []\n"
            "\n\n"
            "calibration.qcal: a histogram over the set of alignment categories\n"
            "(top score, 2nd score, given score)\n"
            "tallying number of alignments that are correct or incorrect.\n"
            "space is currently one of 'G' (genome) or 'T' (transcriptome)\n\n"

            "unscored_fragsort.sam: alignment file sorted by (read id / pair flag) and\n"
            "having alignment score tags given in option -s.\n\n",
            fdef, qdef, Qdef, ldef, Ldef, cdef, sdef, mdef
            );
    return 1;
}


    
struct ScoreCategory
{
    RAW_SCORE_T top_score;
    RAW_SCORE_T sec_score;
    RAW_SCORE_T given_score;

    RawScoreBetter score_better;

    ScoreCategory(RAW_SCORE_T _t, RAW_SCORE_T _s, RAW_SCORE_T _g,
                  bool _lsc) : 
        top_score(_t), sec_score(_s), given_score(_g),
        score_better(_lsc) { }

    ScoreCategory() : 
        top_score(0),
        sec_score(0),
        given_score(0),
        score_better(true) { }

    bool operator<(ScoreCategory const& b) const
    {
        return this->score_better(b.top_score, this->top_score)
            || (this->top_score == b.top_score
                && (this->score_better(b.sec_score, this->sec_score)
                    || (this->sec_score == b.sec_score
                        && this->score_better(b.given_score, this->given_score))));
    }
};


typedef std::map<ScoreCategory, std::pair<RAW_SCORE_T, size_t> > RAW_HISTO;

int main_score_dist(int argc, char ** argv)
{
    char c;
    size_t fdef = 100000;
    float qdef = 0.05;
    float Qdef = 0.95;
    size_t ldef = 0;
    size_t Ldef = 1000000;
    float cdef = 1.0;
    char const* sdef = "AS";
    size_t mdef = 0;

    size_t num_top_fragments_used = fdef;
    float frag_dist_low_quantile = qdef;
    float frag_dist_high_quantile = Qdef;
    size_t min_allowed_fragment_length = ldef;
    size_t max_allowed_fragment_length = Ldef;
    size_t min_fraction_correct_bases = cdef;
    char const* raw_score_tag = sdef;
    bool larger_score_better = false;
    size_t max_valid_fragment_score = mdef;
    bool equivalency_correctness = false;
    char const* expected_read_layout = "";

    while ((c = getopt(argc, argv, "f:q:Q:l:L:c:es:im:y:")) >= 0)
    {
        switch(c)
        {
        case 'f': num_top_fragments_used = static_cast<size_t>(atoi(optarg)); break;
        case 'q': frag_dist_low_quantile = atof(optarg); break;
        case 'Q': frag_dist_high_quantile = atof(optarg); break;
        case 'l': min_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'L': max_allowed_fragment_length = static_cast<size_t>(atoi(optarg)); break;
        case 'c': min_fraction_correct_bases = atof(optarg); break;
        case 'e': equivalency_correctness = true; break;
        case 's': raw_score_tag = optarg; break;
        case 'i': larger_score_better = true; break;
        case 'm': max_valid_fragment_score = static_cast<size_t>(atoi(optarg)); break;
        case 'y': expected_read_layout = optarg; break;

        default: return score_dist_usage(fdef, qdef, Qdef, ldef, Ldef, cdef, sdef, mdef); break;
        }
    }

    if (argc != optind + 2)
    {
        return score_dist_usage(fdef, qdef, Qdef, ldef, Ldef, cdef, sdef, mdef);
    }

    char * unscored_sam_file = argv[optind];
    char * score_calibration_file = argv[optind + 1];

    FILE * unscored_sam_fh = open_or_die(unscored_sam_file, "r", "Input unscored sam file");
    FILE * score_calibration_fh = open_or_die(score_calibration_file, "w", "Output score calibration file");

    SamOrder sam_order(SAM_RID_POSITION, "NONE");

    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(unscored_sam_fh);
    sam_order.AddHeaderContigStats(unscored_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt, expected_read_layout);
    
    //alignments from different physical fragments
    bool paired_reads_are_same_stranded = false;
    bool allow_absent_seq_qual = true; // why not?

    //sources, and take the unique subset of them

    RawScoreBetter score_order(larger_score_better);

    //tally / estimate fragment lengths from top-scoring fragment alignments
    SamBuffer tally_buffer(&sam_order, paired_reads_are_same_stranded);

    size_t min_fragment_length;
    size_t max_fragment_length;

    assert(false);
    //This is currently broken.
    // QuantileFragmentEstimate(min_allowed_fragment_length,
    //                          max_allowed_fragment_length,
    //                          frag_dist_low_quantile,
    //                          frag_dist_high_quantile,
    //                          &unscored_sam_fh, 
    //                          score_order, 
    //                          tally_buffer, 
    //                          raw_score_tag,
    //                          max_valid_fragment_score,
    //                          num_top_fragments_used,
    //                          &min_fragment_length,
    //                          &max_fragment_length);

    bool new_fragment;

    SamBuffer sam_buffer(&sam_order, paired_reads_are_same_stranded);

    char prev_qname[1024] = "";
    std::vector<RAW_SCORE_T> packed_scores;

    RAW_HISTO score_cat_histo;
    read_coords first_guide_coords;
    read_coords second_guide_coords;
    bool seen_a_read = false;
    size_t prev_qid = 0;

    RAW_SCORE_T top_raw_score;
    RAW_SCORE_T sec_raw_score;

    bool guide_is_ones_based_pos = true;

    SamLine * low_bound = NULL;

    while (! feof(unscored_sam_fh))
    {
        NextLine(unscored_sam_fh, sam_buffer, allow_absent_seq_qual,
                 &new_fragment, &seen_a_read, prev_qname, &prev_qid, &low_bound);

        if (new_fragment)
        {
            packed_scores = 
                GetPackedScores(sam_buffer, 
                                min_fragment_length,
                                max_fragment_length,
                                raw_score_tag, 
                                max_valid_fragment_score);


            //find top two scores
            
            std::vector<RAW_SCORE_T> fscopy(packed_scores.size());

            for (size_t p = 0; p != packed_scores.size(); ++p)
            {
                fscopy[p] = UnpackScore(packed_scores[p], max_valid_fragment_score);
            }

            std::sort(fscopy.begin(), fscopy.end(), score_order);
            std::vector<RAW_SCORE_T>::iterator new_end = 
                std::unique(fscopy.begin(), fscopy.end());

            size_t num_unique = std::distance(fscopy.begin(), new_end);
            if (num_unique >= 2)
            {
                top_raw_score = fscopy[0];
                sec_raw_score = fscopy[1];
            }
            else if (num_unique == 1)
            {
                top_raw_score = fscopy[0];
                sec_raw_score = max_valid_fragment_score + 1;
            }
            else
            {
                //this is really unnecessary since it means the
                //sam_buffer doesn't have any entries.
                top_raw_score = max_valid_fragment_score + 1;
                sec_raw_score = max_valid_fragment_score + 1;
            }

            //then, for each fragment, evaluate whether it is correctly aligned,
            //and update the histogram

            PAIRED_READ_SET::const_iterator pit = sam_buffer.unique_entry_pairs.begin();

            size_t frag_index = 0;
            for (pit = sam_buffer.unique_entry_pairs.begin();
                 pit != sam_buffer.unique_entry_pairs.end(); ++pit)
            {
                SamLine const* first = (*pit).first;
                SamLine const* second = (*pit).second;

                CigarFromSimSAMLine(first->qname, 
                                    first->first_read_in_pair(),
                                    guide_is_ones_based_pos,
                                    &first_guide_coords);
                
                CigarFromSimSAMLine(second->qname, 
                                    second->first_read_in_pair(),
                                    guide_is_ones_based_pos,
                                    &second_guide_coords);
                
                //calculate whether the alignment is as good as the correct one.
                RAW_SCORE_T fragment_score = 
                    UnpackScore(packed_scores[frag_index], max_valid_fragment_score);

                ScoreCategory score_cat(top_raw_score, sec_raw_score,
                                        fragment_score,
                                        larger_score_better);

                if (score_cat_histo.find(score_cat) == score_cat_histo.end())
                {
                    score_cat_histo[score_cat] = std::make_pair(0, 0);
                }

                bool is_correct;

                if (equivalency_correctness)
                {
                    //consider fragment alignments 'correct' if the
                    //aligner found the same or better raw alignment
                    //score as the correct one.
                    RAW_SCORE_T correct_fragment_score = 
                    first_guide_coords.num_errors +
                    second_guide_coords.num_errors;
                    
                    //correct if the correct fragment score is not better than the given one.
                    is_correct = 
                        (! score_order(correct_fragment_score, fragment_score));
                }
                else
                {
                    //use the traditional notion of alignment
                    //correctness as that the aligner correctly placed
                    //enough bases
                    Cigar::CIGAR_VEC first_guide_cigar = 
                        Cigar::FromString(first_guide_coords.cigar, first_guide_coords.position);
                
                    std::multiset<std::pair<size_t, size_t> > first_guide_cigar_index =
                        Cigar::ComputeOffsets(first_guide_cigar);
                
                    Cigar::CIGAR_VEC second_guide_cigar = 
                        Cigar::FromString(second_guide_coords.cigar, second_guide_coords.position);
                
                    std::multiset<std::pair<size_t, size_t> > second_guide_cigar_index =
                        Cigar::ComputeOffsets(second_guide_cigar);

                    size_t first_num_guide, first_num_test;
                    size_t first_num_correct_bases = 
                        CountCorrectBases(first, first_guide_coords,
                                          first_guide_cigar,
                                          first_guide_cigar_index,
                                          &first_num_guide, &first_num_test);
                
                    size_t second_num_guide, second_num_test;
                    size_t second_num_correct_bases = 
                        CountCorrectBases(second, second_guide_coords,
                                          second_guide_cigar,
                                          second_guide_cigar_index,
                                          &second_num_guide, &second_num_test);
                
                
                    size_t num_correct_bases = first_num_correct_bases + second_num_correct_bases;
                    size_t num_guide_bases = first_num_guide + second_num_guide;

                    float frac_correct_bases = num_correct_bases > 0
                        ? static_cast<float>(num_correct_bases) / static_cast<float>(num_guide_bases)
                        : 0.0;

                    is_correct = frac_correct_bases >= min_fraction_correct_bases;
                }

                if (is_correct)
                {
                    score_cat_histo[score_cat].first++;
                }
                else
                {
                    score_cat_histo[score_cat].second++;
                }

                ++frag_index;
            }

            sam_buffer.purge(NULL, NULL, NULL, low_bound);
        }
    }

    RAW_HISTO::const_iterator rit;
    fprintf(score_calibration_fh, "score_tag: %s\n", raw_score_tag);

    for (rit = score_cat_histo.begin(); rit != score_cat_histo.end(); ++rit)
    {
        ScoreCategory const& cat = (*rit).first;
        size_t num_right = (*rit).second.first;
        size_t num_wrong = (*rit).second.second;

        fprintf(score_calibration_fh, "%i\t%i\t%i\t%Zu\t%Zu\n",
                cat.top_score,
                cat.sec_score,
                cat.given_score,
                num_right, num_wrong);
    }

    fclose(score_calibration_fh);
    fclose(unscored_sam_fh);

    return 0;
}
