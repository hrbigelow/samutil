#ifndef _ALIGN_EVAL_AUX_H
#define _ALIGN_EVAL_AUX_H

#include <vector>

#include "sam_helper.h"
#include "sam_order.h"

struct LineIndex
{
    size_t index;
    size_t start_offset;
    int line_length;
    LineIndex(size_t _i, size_t _s, int _l) : index(_i), start_offset(_s), line_length(_l) { }
    LineIndex() : index(0), start_offset(0), line_length(0) { }
};


typedef std::vector<LineIndex>::iterator INDEX_ITER;

bool less_offset(LineIndex const& a, LineIndex const& b);
bool less_key(LineIndex const& a, LineIndex const& b);


struct partial_index_aux
{
    SamOrder const* sam_order;
    partial_index_aux(SamOrder const* _sam_order);
    LineIndex operator()(char * samline);
};


// input: full_samline raw string
// output: in-place modified string containing truncated string
struct truncate_sam_unary
{
    truncate_sam_unary();
    char * operator()(char * full_samline);
};


std::pair<size_t, size_t> 
process_chunk(char * chunk_buffer_in,
              char * chunk_buffer_out,
              size_t chunk_length,
              SamOrder const& sam_order,
              FILE * chunk_tmp_fh,
              std::vector<LineIndex> * line_index);

void
get_key_quantiles(std::vector<LineIndex> const& line_index,
                  size_t num_chunks,
                  size_t * key_quantile_sizes,
                  std::vector<LineIndex> * key_quantile_sentinels);


void 
write_final_merge(std::vector<LineIndex> const& ok_index,
                  std::vector<INDEX_ITER> const& offset_quantiles,
                  FILE * tmp_fhs[],
                  size_t num_chunks,
                  FILE * out_dat_fh,
                  FILE * out_ind_fh);


std::vector<INDEX_ITER> 
get_quantiles(std::vector<LineIndex> * line_index,
              bool (less_fcn)(LineIndex const&, LineIndex const&),
              size_t num_chunks);


#endif // _ALIGN_EVAL_AUX_H
