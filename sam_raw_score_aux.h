#ifndef _SAM_RAW_SCORE_AUX_H
#define _SAM_RAW_SCORE_AUX_H

#include "sam_buffer.h"
#include "align_eval_raw.h"
#include "cigar_ops.h"

struct ScorePair
{
    size_t top_score;
    size_t sec_score;
    ScorePair(size_t _t, size_t _s);
    ScorePair();
    bool operator<(ScorePair const& s) const;
};


size_t FragmentScore(size_t top_raw_score, 
                     size_t sec_raw_score,
                     size_t default_score_if_missing,
                     int template_length,
                     size_t min_allowed_template_length,
                     size_t max_allowed_template_length);

ScorePair FindTopTwoScores(SamBuffer const& sam_buffer, 
                           size_t min_allowed_fragment_length,
                           size_t max_allowed_fragment_length,
                           char const* tag,
                           size_t missing_score_default,
                           bool larger_is_better,
                           std::vector<size_t> * fragment_scores);


size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index);


void NextLine(FILE * unnscored_sam_fh, 
              SamBuffer & sam_buffer,
              bool ones_based_pos,
              bool allow_absent_seq_qual,
              bool * new_fragment, 
              bool * ignore_sambuffer_bound,
              bool * seen_a_read,
              char * prev_rname,
              size_t * prev_qid);


#endif // _SAM_RAW_SCORE_AUX_H
