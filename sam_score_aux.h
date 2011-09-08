#ifndef _SAM_RAW_SCORE_AUX_H
#define _SAM_RAW_SCORE_AUX_H

#include <map>
#include <utility>

#include "sam_buffer.h"
#include "sam_helper.h"
#include "align_eval_raw.h"
#include "cigar_ops.h"

typedef int RAW_SCORE_T;



struct ScoreSpace
{
    int raw_score;
    char alignment_space;
public:
    ScoreSpace(int _rs, char _as) : raw_score(_rs), alignment_space(_as) { }
};


struct RawScoreComp
{
    bool larger_score_better;
    RawScoreComp(bool _lsb = true) : larger_score_better(_lsb) { }

    //returns true if a is worse than b
    bool operator()(int a, int b) const
    {
        return this->larger_score_better ? a < b : a > b;
    }
};



class FragmentScore
{

 public:
    size_t min_fragment_length;
    size_t max_fragment_length;
    int min_fragment_score;
    int max_fragment_score;
    int worst_fragment_score; //min or max depending on logic
    bool larger_score_better;
    char raw_score_tag[10];
    char * space_ordering; // ordering from best to worst
    CONTIG_SPACE contig_space; //map def
    int * score_table; //
    size_t table_shift_bits; // for use with score_table_index

    RawScoreComp raw_score_comp;

    FragmentScore(size_t _minfrag, size_t _maxfrag) :
        min_fragment_length(_minfrag),
        max_fragment_length(_maxfrag),
        space_ordering(NULL),
        score_table(NULL) { }

    ~FragmentScore();

    void init(char const* qcal_file, char const* contig_space_file);

    int raw_score(SamLine const* samline) const;
    char alignment_space(SamLine const* samline) const;
    size_t score_table_index(int top_score, int sec_score, int given_score) const;
    bool operator()(ScoreSpace const& a, ScoreSpace const& b) const;
};


class FragmentScoreWrap
{
    FragmentScore const* fragment_score;
public:
    FragmentScoreWrap(FragmentScore const* fs);
    bool operator()(ScoreSpace const& a, ScoreSpace const& b) const;

};


void set_score_fields(SamBuffer const& sam_buffer, 
                      FragmentScore const& fragment_scoring);


size_t CountCorrectBases(SamLine const* samline, 
                         read_coords const& guide_coords, 
                         Cigar::CIGAR_VEC const& guide_cigar, 
                         Cigar::CIGAR_INDEX const& guide_cigar_index,
                         size_t * num_guide_bases,
                         size_t * num_test_bases);


void NextLine(FILE * unnscored_sam_fh, 
              SamBuffer & sam_buffer,
              bool * new_fragment, 
              bool * seen_a_read,
              char * prev_rname,
              size_t * prev_qid,
              SamLine ** low_bound);


/*
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
*/


#endif // _SAM_RAW_SCORE_AUX_H
