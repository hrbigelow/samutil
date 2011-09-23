#ifndef _SAM_AUX_H
#define _SAM_AUX_H

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include "sam_score_aux.h"
#include "align_eval_aux.h"


typedef std::set<SequenceProjection>::const_iterator SP_ITER;

#ifdef __GXX_EXPERIMENTAL_CXX0X__
typedef std::unordered_map<char const*, SP_ITER, to_integer, eqstr> PROJ_MAP;
#else
typedef std::tr1::unordered_map<char const*, SP_ITER, to_integer, eqstr> PROJ_MAP;
#endif

typedef PROJ_MAP::const_iterator NP_ITER;


bool tx_nonoverlapping_or_passthrough(PROJ_MAP const& tx_projections,
                                      char const* tx1, char const* tx2);


// inputs: SamLine * [beg, end) range, SamBuffer *
// outputs: char * projected_deduped_samlines
// allocates output.
struct project_dedup_print
{
    PROJ_MAP const* projections;
    SamOrder const* tx_sam_order;
    SamOrder const* genome_sam_order;
    bool inserts_are_introns;
    size_t average_rsam_line_length;

    project_dedup_print(PROJ_MAP const* _proj_map,
                        SamOrder const* _tx_order,
                        SamOrder const* _genome_order,
                        bool _inserts_are_introns,
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
