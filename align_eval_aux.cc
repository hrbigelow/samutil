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
    LineIndex operator()(char const* samline)
    {
        return LineIndex((this->sam_order->*(this->sam_order->sam_index))(samline), 
                         0,
                         strlen(samline) + 1);
    }        
};



std::vector<LineIndex> 
build_index(char const* sam_file,
            char * chunk_buffer,
            size_t max_mem,
            size_t max_line,
            SamOrder const& sam_order,
            // size_t (* samline_pos)(char const*), 
            size_t * num_chunks)
{

    std::vector<LineIndex> line_index;

    FILE * sam_fh = fopen(sam_file, "r");
    SetToFirstDataLine(&sam_fh);
    size_t header_length = ftell(sam_fh);
    
    gzFile sam_zfh = gzopen(sam_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(sam_zfh, max_mem, max_line);
    gzclose(sam_zfh);

    chunk_lengths[0] -= header_length;
    *num_chunks = chunk_lengths.size();

    // size_t max_chunk_length = 
    //     std::max_element(chunk_lengths.begin(), chunk_lengths.end());

    size_t current_offset = header_length;
    fprintf(stderr, "chunks [0-%zu] ", *num_chunks - 1);
    fflush(stderr);

    for (size_t chunk = 0; chunk != *num_chunks; ++chunk)
    {
        size_t nbytes_read = fread(chunk_buffer, 1, chunk_lengths[chunk], sam_fh);
        assert(nbytes_read == chunk_lengths[chunk]);
        chunk_buffer[nbytes_read] = '\0';

        //this is really fast
        std::vector<char *> samlines = FileUtils::find_line_starts(chunk_buffer);

        size_t S = samlines.size();

        LineIndex * line_index_chunk = new LineIndex[S];

        //this is really fast
        char * line_end = chunk_buffer;
        while ((line_end = strchr(line_end, '\n')) != NULL)
        {
            *line_end++ = '\0';
        }

        std::vector<char *>::iterator sit;
        LineIndex * lit;
        
        partial_index_aux paux(&sam_order);

        __gnu_parallel::transform(samlines.begin(), samlines.end(), line_index_chunk, paux);

        // now update the current offsets
        for (lit = line_index_chunk; lit != line_index_chunk + S; ++lit)
        {
            (*lit).start_offset = current_offset;
            current_offset += (*lit).line_length;
        }

        line_index.insert(line_index.end(), line_index_chunk, line_index_chunk + S);

        delete line_index_chunk;

        fprintf(stderr, "%zu ", chunk);
        fflush(stderr);
    }
    fclose(sam_fh);

    return line_index;
}
