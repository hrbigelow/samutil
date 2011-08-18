#ifndef _SAM_RAW_SCORE_AUX_H
#define _SAM_RAW_SCORE_AUX_H

#include <map>
#include <utility>

#include "sam_buffer.h"
#include "align_eval_raw.h"
#include "cigar_ops.h"

typedef int RAW_SCORE_T;

struct RawScoreBetter
{
    bool larger_score_better;
    RawScoreBetter(bool _lsb) : larger_score_better(_lsb) { }

    //returns true if a is 'better' than b
    bool operator()(RAW_SCORE_T a, RAW_SCORE_T b) const
    {
        return this->larger_score_better ? a > b : a < b;
    }
};


union RawScoreTriplet
{
    size_t data;
    struct
    {
        unsigned int top_raw_score : 16;
        unsigned int sec_raw_score : 16;
        unsigned int given_raw_score : 16;
        int pad : 16;
    } f;
};


struct CollapseScoreTriplet
{
    size_t min_score;
    size_t max_score;
    size_t shift_bits;

    CollapseScoreTriplet(size_t _min, size_t _max);
    size_t operator()(size_t top, size_t sec, size_t given) const;
};



std::vector<RAW_SCORE_T> 
GetPackedScores(SamBuffer const& sam_buffer, 
                size_t min_allowed_fragment_length,
                size_t max_allowed_fragment_length,
                char const* raw_score_tag,
                RAW_SCORE_T missing_default_score);


RAW_SCORE_T UnpackScore(RAW_SCORE_T packed_score,
                        RAW_SCORE_T max_valid_fragment_score);

void SetScoreFields(SamBuffer const& sam_buffer, 
                    size_t min_allowed_fragment_length,
                    size_t max_allowed_fragment_length,
                    char const* raw_score_tag,
                    RAW_SCORE_T missing_default_score,
                    bool larger_raw_score_is_better,
                    int const* score_cpd,
                    CollapseScoreTriplet const& score_triplet);


size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index,
                         size_t * num_guide_bases,
                         size_t * num_test_bases);


void NextLine(FILE * unnscored_sam_fh, 
              SamBuffer & sam_buffer,
              bool allow_absent_seq_qual,
              bool * new_fragment, 
              bool * seen_a_read,
              char * prev_rname,
              size_t * prev_qid,
              SamLine ** low_bound);


std::vector<size_t>
TallyFragmentLengths(FILE ** sam_fh, 
                     RawScoreBetter const& score_order,
                     SamBuffer & tally_buffer,
                     char const* score_tag,
                     RAW_SCORE_T missing_default_score,
                     size_t num_top_fragments_used);


//estimate fragment lengths based on quantiles
void QuantileFragmentEstimate(size_t min_allowed_fragment_length,
                              size_t max_allowed_fragment_length,
                              float frag_dist_low_quantile,
                              float frag_dist_high_quantile,
                              FILE ** sam_fh,
                              RawScoreBetter const& score_order,
                              SamBuffer & tally_buffer,
                              char const* score_tag,
                              RAW_SCORE_T missing_default_score,
                              size_t num_top_scoring_frags_used,
                              size_t * min_est_fragment_length, 
                              size_t * max_est_fragment_length);


CollapseScoreTriplet
ParseScoreCalibration(FILE * score_cal_fh,
                      char * raw_score_tag,
                      bool * larger_score_is_better,
                      int * & score_cpd);


#endif // _SAM_RAW_SCORE_AUX_H
