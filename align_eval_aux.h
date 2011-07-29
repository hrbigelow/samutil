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

std::vector<LineIndex> 
build_index(char const* sam_file,
            char * chunk_buffer,
            size_t max_mem,
            size_t max_line,
            SamOrder const& sam_order,
            /* size_t (* samline_pos)(char const*),  */
            size_t * num_chunks);


#endif // _ALIGN_EVAL_AUX_H
