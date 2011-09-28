#ifndef _SAM_AUX_H
#define _SAM_AUX_H

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include "sam_score_aux.h"
#include "align_eval_aux.h"


typedef std::set<SequenceProjection, less_seq_projection>::const_iterator SEQ_PROJ_ITER;

#ifdef __GXX_EXPERIMENTAL_CXX0X__
typedef std::unordered_map<char const*, SEQ_PROJ_ITER, to_integer, eqstr> PROJ_MAP;
#else
typedef std::tr1::unordered_map<char const*, SEQ_PROJ_ITER, to_integer, eqstr> PROJ_MAP;
#endif

typedef PROJ_MAP::const_iterator NP_ITER;


SEQ_PROJ_ITER 
non_overlapping_range(SEQ_PROJ_ITER start, SEQ_PROJ_ITER set_end);


// inputs: SamLine * [beg, end) range, SamBuffer *
// outputs: char * projected_deduped_samlines
// allocates output.
struct project_dedup_print
{
    PROJ_MAP const* projections;
    SamOrder const* tx_sam_order;
    SamOrder const* genome_sam_order;
    bool inserts_are_introns;
    bool retain_unsequenced_projection;
    size_t average_rsam_line_length;

    project_dedup_print(PROJ_MAP const* _proj_map,
                        SamOrder const* _tx_order,
                        SamOrder const* _genome_order,
                        bool _inserts_are_introns,
                        bool _retain_unseq,
                        size_t _ave_len);

    std::vector<char> * operator()(std::pair<SAMIT, SAMIT> const& range);

};


// inputs: rSAM-formatted SamLine *, INDEX_ITER
// outputs: SAM-formatted allocated SAM record string
struct rsam_to_sam_binary
{
    char const* seq_buffer; // generated from 'samutil seqindex'
    size_t data_buffer_offset; // the start_offset corresponding to the beginning of this seq_buffer

    rsam_to_sam_binary(char const* _seq_buffer, size_t _dbo);
    char * operator()(SamLine * samline, INDEX_ITER li_iter);

};

#endif // _SAM_AUX_H
