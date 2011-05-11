#ifndef _SAM_RAW_SCORE_AUX_H
#define _SAM_RAW_SCORE_AUX_H

#include "sam_buffer.h"
#include "align_eval_raw.h"
#include "cigar_ops.h"

enum AlignmentSpace
    {
        GENOME,
        TRANSCRIPTOME
    };


#define AlignmentSpaces "GT"
#define AlignmentSpaceMissingDefault 'G'


AlignmentSpace ParseAlignmentSpace(char as);
char PrintAlignmentSpace(AlignmentSpace const& as);


struct SpaceScore
{
    size_t score;
    AlignmentSpace space;
    SpaceScore(size_t _s, AlignmentSpace _as);
    SpaceScore();
    bool operator==(SpaceScore const& b) const;
    bool operator<(SpaceScore const& b) const;
};



struct SpaceScore_better
{
    bool larger_score_better;
    bool larger_space_better;
    bool operator()(SpaceScore const& a, SpaceScore const& b) const;
    SpaceScore_better(bool _lscore, bool _lspace);
};



//holds the pair of top and second score, used for classifying
//the entire set of alignments of a given fragment.
struct ScorePair
{
    SpaceScore top_score;
    SpaceScore sec_score;
    ScorePair(SpaceScore _t, SpaceScore _s);
    ScorePair();

    //needed just for the SCORE_CPD object.
    bool operator<(ScorePair const& s) const;
};


typedef std::vector<SpaceScore> ScoreVec;

SpaceScore FragmentScore(size_t top_raw_score, 
                         size_t sec_raw_score,
                         AlignmentSpace space,
                         size_t missing_default_score,
                         size_t template_length,
                         size_t min_allowed_template_length,
                         size_t max_allowed_template_length);

ScorePair FindTopTwoScores(SamBuffer const& sam_buffer, 
                           size_t min_allowed_fragment_length,
                           size_t max_allowed_fragment_length,
                           char const* tag,
                           size_t missing_default_score,
                           SpaceScore_better const& score_better_cmp,
                           ScoreVec * fragment_scores);


size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index);


void NextLine(FILE * unnscored_sam_fh, 
              SamBuffer & sam_buffer,
              bool allow_absent_seq_qual,
              bool * new_fragment, 
              bool * ignore_sambuffer_bound,
              bool * seen_a_read,
              char * prev_rname,
              size_t * prev_qid);


std::vector<double>
TallyFragmentLengths(FILE ** sam_fh, 
                     SpaceScore_better const& score_best_is_less,
                     SamBuffer & tally_buffer,
                     char const* score_tag,
                     size_t missing_default_score,
                     float fraction_top_scoring_frags_used);

#endif // _SAM_RAW_SCORE_AUX_H
