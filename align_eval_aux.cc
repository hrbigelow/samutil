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
// chunk_buffer_in is set to 
//line_index will be appended with a K-sorted chunk range
/*
1. index
2. sort
3. write sorted chunk to buffer
4. write buffer to file
5. append sorted index to main index

Returns a pair of <number bytes, number of lines> in chunk processed
*/

    // size_t nbytes_read = fread(chunk_buffer_in, 1, chunk_length, input_sam_fh);
    // assert(nbytes_read == chunk_length);
    // chunk_buffer_in[nbytes_read] = '\0';

std::pair<size_t, size_t> 
process_chunk(char * chunk_buffer_in,
              char * chunk_buffer_out,
              size_t chunk_length,
              SamOrder const& sam_order,
              FILE * chunk_tmp_fh,
              std::vector<LineIndex> * line_index)
{

    //this is really fast
    std::vector<char *> samlines = FileUtils::find_line_starts(chunk_buffer_in);

    size_t S = samlines.size();

    if (S == 0)
    {
        return std::pair<size_t, size_t>(0, 0);
    }

    LineIndex * line_index_chunk = new LineIndex[S];

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

    size_t nbytes_written = std::distance(chunk_buffer_out, write_pointer);
    assert(nbytes_written == chunk_length);

    fwrite(chunk_buffer_out, 1, chunk_length, chunk_tmp_fh);
    fflush(chunk_tmp_fh);

    (*line_index).insert((*line_index).end(), line_index_chunk, line_index_chunk + S);
    delete line_index_chunk;

    return std::make_pair(nbytes_written, S);
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


// gathers next sub-chunks from each temp file, starting at
// 'prev_chunk_starts', and using 'query_key_quantile' to find
// the end of each subchunk.
// postconditions: 
// catenated_index populated and ordered as (O, K)
// chunk_buffer populated to match catenated_index
// prev_chunk_starts updated to next set of chunk starts for next call
// returns number of bytes written
size_t catenate_subchunks(LineIndex const& query_key_quantile,
                          FILE * tmp_fhs[],
                          std::vector<INDEX_ITER> * prev_chunk_starts,
                          std::vector<INDEX_ITER> const& chunk_ends, /* ends of each chunk ordered by (O, K) */
                          size_t num_chunks,
                          char ** chunk_buffer,
                          std::vector<LineIndex> * catenated_index,
                          size_t ** subrange_sizes)
{
    std::pair<INDEX_ITER, INDEX_ITER> * subrange_iters =
        new std::pair<INDEX_ITER, INDEX_ITER>[num_chunks];

    char * write_pointer = (*chunk_buffer);
    size_t buffer_pos = 0;
    size_t total_subrange_size = 0;

    INDEX_ITER beg, end, sub_k;
    size_t S;

    for (size_t o = 0; o != num_chunks; ++o)
    {
        beg = (*prev_chunk_starts)[o];
        end = chunk_ends[o+1];
        sub_k = std::lower_bound(beg, end, query_key_quantile, less_key);

        //assert(__gnu_cxx::is_sorted(beg, end, less_key));

        // fprintf(stdout, "key %Zu, offset %Zu, num_lines: %Zu\n",
        //         k, o, part);

        S = 0;
        for (INDEX_ITER kit = beg; kit != sub_k; ++kit)
        {
            (*kit).start_offset = S + buffer_pos;
            S += (*kit).line_length;
        }
        fread(write_pointer, 1, S, tmp_fhs[o]);
        write_pointer += S;
        buffer_pos += S;

        // last_merged_size = merged.size();
        subrange_iters[o] = std::make_pair(beg, sub_k);
        (*subrange_sizes)[o] = std::distance(beg, sub_k);
        total_subrange_size += (*subrange_sizes)[o];

        // fprintf(stdout, "key %Zu, offset %Zu, added: %Zu, num_lines: %Zu\n",
        //         k, o, part, merged.size());

        (*prev_chunk_starts)[o] = sub_k;
    }

    (*catenated_index).reserve(total_subrange_size);
    (*catenated_index).resize(0);

    for (size_t o = 0; o != num_chunks; ++o)
    {
        catenated_index->insert(catenated_index->end(), 
                                subrange_iters[o].first, 
                                subrange_iters[o].second);
    }

    delete subrange_iters;
    return buffer_pos; // number of bytes written
}


// merge an (O, K) ordered index.
// precondition: ok_index is an (O, K) ordered index in 'num_ranges' pieces,
// of size 'subrange_sizes'
// postcondition: k_index is the (K) ordered index
// warning: destroys the original ok_index

void merge_ok_index(size_t const* subrange_sizes,
                    size_t num_ranges,
                    std::vector<LineIndex> ** ok_index,
                    std::vector<LineIndex> ** k_index)
{

    size_t * sz_tmp = new size_t[num_ranges + 1];
    std::copy(subrange_sizes, subrange_sizes + num_ranges, sz_tmp);
    sz_tmp[num_ranges] = 0;

    size_t total_size = std::accumulate(subrange_sizes, subrange_sizes + num_ranges, 0);

    (*k_index)->reserve(total_size);
    (*k_index)->resize(total_size);

    //here's the tricky part.
    //iteratively merge consecutive ranges until there is only one.
    //invariant: subrange_sizes valid, num_ranges valid, k_index valid
    size_t num_merge_ranges = num_ranges + (num_ranges % 2); // square it off

    std::vector<LineIndex>::iterator swapit, mergit;

    while (num_merge_ranges > 1)
    {
        std::swap((*k_index), (*ok_index));
        swapit = (*k_index)->begin();
        mergit = (*ok_index)->begin();

        for (size_t o = 0; o != num_merge_ranges; o += 2)
        {
            __gnu_parallel::merge(swapit, swapit + sz_tmp[o],
                                  swapit + sz_tmp[o],
                                  swapit + sz_tmp[o] + sz_tmp[o + 1],
                                  mergit, less_key);
                    
            swapit += sz_tmp[o] + sz_tmp[o + 1];
            mergit += sz_tmp[o] + sz_tmp[o + 1];

            sz_tmp[o / 2] = sz_tmp[o] + sz_tmp[o + 1];
        }
        num_merge_ranges /= 2;

        if (num_merge_ranges > 1 && (num_merge_ranges % 2) != 0)
        {
            num_merge_ranges += (num_merge_ranges % 2);
            sz_tmp[num_merge_ranges - 1] = 0;
        }
    }
    delete sz_tmp;

}

//write blocks from 'unordered' to 'ordered' according to 'ordering'
//index.  return number of bytes written
size_t write_new_ordering(char const* unordered, 
                          std::vector<LineIndex> * ordering,
                          char * ordered)
{
    char * write_pointer = ordered;
    
    std::vector<LineIndex>::iterator mit;
    size_t S = 0;

    for (mit = (*ordering).begin(); mit != (*ordering).end(); ++mit)
    {
        memcpy(write_pointer, unordered + (*mit).start_offset,
               (*mit).line_length);
        
        write_pointer += (*mit).line_length;
        S += (*mit).line_length;
        (*mit).start_offset = S;
    }
    return std::distance(ordered, write_pointer);
}


size_t set_start_offsets(std::vector<LineIndex>::iterator beg,
                         std::vector<LineIndex>::iterator end,
                         size_t initial_offset)
{
    std::vector<LineIndex>::iterator cur;
    size_t off = initial_offset;
    for (cur = beg; cur != end; ++cur)
    {
        (*cur).start_offset = off;
        off += (*cur).line_length;
    }
    return off;
}
                         


// Preconditions: 
// ok_index: sorted by (O, K)
// tmp_fhs: opened and at beginning. hold K-sorted chunks consistent with ok_index
// offset_quantiles define positions in ok_index corresponding to tmp_fhs chunks
// Postconditions:
// out_dat_fh has final K-sorted records
// out_ind_fh has final K-sorted index entries

void write_final_merge(std::vector<LineIndex> const& ok_index,
                       std::vector<INDEX_ITER> const& offset_quantiles,
                       FILE * tmp_fhs[],
                       size_t num_chunks,
                       FILE * out_dat_fh,
                       FILE * out_ind_fh)
{
    std::vector<LineIndex> key_quantile_sentinels;
    size_t * key_quantile_sizes = new size_t[num_chunks];

    // 2. Find N quantile key values in index.  does not change
    // order of line_index.
    get_key_quantiles(ok_index, num_chunks, key_quantile_sizes, & key_quantile_sentinels);

    size_t kq_size_max = 
        *std::max_element(key_quantile_sizes, key_quantile_sizes + num_chunks);

    // 4. Allocate chunk of memory to hold biggest key_quantile
    char * chunk_buffer_in = new char[kq_size_max];
    char * chunk_buffer_out = new char[kq_size_max];

    size_t * subrange_sizes = new size_t[num_chunks];

    std::vector<INDEX_ITER> prev_sub_k(offset_quantiles.size() - 1);
    std::copy(offset_quantiles.begin(), offset_quantiles.end() - 1, prev_sub_k.begin());

    fprintf(stderr, "Merging chunks [0-%zu]:", num_chunks - 1);
    fflush(stderr);

    size_t chunk_offset = 0;

    for (size_t k = 0; k != num_chunks; ++k)
    {

        std::vector<LineIndex> buf1, buf2;
        std::vector<LineIndex> * ok_index_ptr = & buf1;
        std::vector<LineIndex> * k_index_ptr = & buf2;

        size_t bytes_written = 
            catenate_subchunks(key_quantile_sentinels[k],
                               tmp_fhs,
                               & prev_sub_k,
                               offset_quantiles,
                               num_chunks,
                               & chunk_buffer_in,
                               ok_index_ptr,
                               & subrange_sizes);
        
        assert(bytes_written == key_quantile_sizes[k]);

        merge_ok_index(subrange_sizes, num_chunks, & ok_index_ptr, & k_index_ptr);

        bytes_written = 
            write_new_ordering(chunk_buffer_in, k_index_ptr, chunk_buffer_out);

        assert(bytes_written == key_quantile_sizes[k]);

        fwrite(chunk_buffer_out, 1, key_quantile_sizes[k], out_dat_fh);

        if (out_ind_fh != NULL)
        {
            for (std::vector<LineIndex>::const_iterator i = (*k_index_ptr).begin();
                 i != (*k_index_ptr).end(); ++i)
            {
                fprintf(out_ind_fh, "%zu\t%zu\t%i\n",
                        (*i).index, (*i).start_offset + chunk_offset, (*i).line_length);
            }
        }

        chunk_offset += bytes_written;

        fprintf(stderr, " %zu", k);
        fflush(stderr);
    }

    delete chunk_buffer_in;
    delete chunk_buffer_out;
    delete key_quantile_sizes;
    delete subrange_sizes;

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
