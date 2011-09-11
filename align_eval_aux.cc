#include "align_eval_aux.h"
#include "file_utils.h"

#include <parallel/algorithm>

//partial ordering
bool less_offset(LineIndex const& a, LineIndex const& b)
{
    return a.start_offset < b.start_offset;
}


//partial ordering
bool less_key(LineIndex const& a, LineIndex const& b)
{
    return a.index < b.index ||
        (a.index == b.index &&
         a.start_offset < b.start_offset);
}


//unary function for computing index
struct partial_index_aux
{
    SamOrder const* sam_order;
    partial_index_aux(SamOrder const* _sam_order) : sam_order(_sam_order) { }
    LineIndex operator()(char * samline)
    {
        char * end = strchr(samline, '\n');
        *end = '\0';

        LineIndex index((this->sam_order->*(this->sam_order->sam_index))(samline), 
                        0, strlen(samline) + 1);
        *end = '\n';
        return index;
    }        
};


//tmp_fhs are open. 
//chunk_buffer_in and chunk_buffer_out are allocated to
//largest chunk length in chunk_lengths. 
//input_sam_fh is set to first data line
//line_index will be appended with a K-sorted chunk range
/*
1. index
2. sort
3. write sorted chunk to buffer
4. write buffer to file
5. append sorted index to main index

Returns a pair of <number bytes, number of lines> in chunk processed
*/

std::pair<size_t, size_t> 
process_chunk(char * chunk_buffer_in,
              char * chunk_buffer_out,
              size_t chunk_length,
              SamOrder const& sam_order,
              FILE * input_sam_fh,
              FILE * chunk_tmp_fh,
              std::vector<LineIndex> * line_index)
{
    size_t nbytes_read = fread(chunk_buffer_in, 1, chunk_length, input_sam_fh);
    assert(nbytes_read == chunk_length);
    chunk_buffer_in[nbytes_read] = '\0';

    //this is really fast
    std::vector<char *> samlines = FileUtils::find_line_starts(chunk_buffer_in);

    size_t S = samlines.size();

    LineIndex * line_index_chunk = new LineIndex[S];

    //this is really fast
    // char * line_end = chunk_buffer_in;
    // while ((line_end = strchr(line_end, '\n')) != NULL)
    // {
    //     *line_end++ = '\0';
    // }

    std::vector<char *>::iterator sit;
    LineIndex * lit;
        
    partial_index_aux paux(&sam_order);

    __gnu_parallel::transform(samlines.begin(), samlines.end(), line_index_chunk, paux);

    // now update the current offsets
    size_t current_offset = 0;
    for (lit = line_index_chunk; lit != line_index_chunk + S; ++lit)
    {
        (*lit).start_offset = current_offset;
        current_offset += (*lit).line_length;
    }

    __gnu_parallel::sort(line_index_chunk, line_index_chunk + S, less_key);

    char * write_pointer = chunk_buffer_out;
    for (lit = line_index_chunk; lit != line_index_chunk + S; ++lit)
    {
        memcpy(write_pointer,
               chunk_buffer_in + (*lit).start_offset, (*lit).line_length);

        write_pointer += (*lit).line_length;
    }

    assert(static_cast<size_t>(std::distance(chunk_buffer_out, write_pointer)) 
           == chunk_length);

    fwrite(chunk_buffer_out, 1, chunk_length, chunk_tmp_fh);
    fflush(chunk_tmp_fh);

    (*line_index).insert((*line_index).end(), line_index_chunk, line_index_chunk + S);
    delete line_index_chunk;

    return std::make_pair(nbytes_read, S);
}


void
get_key_quantiles(std::vector<LineIndex> const& line_index,
                  size_t num_chunks,
                  size_t * key_quantile_sizes,
                  std::vector<LineIndex> * key_quantile_sentinels)
{
    std::vector<LineIndex> line_index_copy(line_index);

    std::vector<INDEX_ITER> key_quantiles =
        get_quantiles(& line_index_copy, less_key, num_chunks);
    
    std::fill(key_quantile_sizes, key_quantile_sizes + num_chunks, 0);

    for (size_t k = 0; k != num_chunks; ++k)
    {
        INDEX_ITER kit, kit_end = key_quantiles[k+1];

        LineIndex sentinel = (kit_end == line_index_copy.end())
            ? LineIndex(INT64_MAX, INT64_MAX, 0)
            : *kit_end;

        (*key_quantile_sentinels).push_back(sentinel);

        size_t * kq_size = key_quantile_sizes + k;
        for (kit = key_quantiles[k]; kit != kit_end; ++kit)
        {
            (*kq_size) += (*kit).line_length;
        }
    }
}

//what does this do?  

std::vector<INDEX_ITER> 
get_quantiles(std::vector<LineIndex> * line_index,
              bool (less_fcn)(LineIndex const&, LineIndex const&),
              size_t num_chunks)
{
    
    size_t lines_per_chunk = (*line_index).size() / num_chunks;
    std::vector<INDEX_ITER> quantiles;
    
    INDEX_ITER iter = (*line_index).begin();
    INDEX_ITER end = (*line_index).end();

    quantiles.push_back(iter);
    for (size_t n = 0; n != num_chunks - 1; ++n)
    {
        std::advance(iter, lines_per_chunk);
        __gnu_parallel::nth_element(quantiles.back(), iter, end, less_fcn);
        quantiles.push_back(iter);
    }
    quantiles.push_back(end);
    return quantiles;
}
