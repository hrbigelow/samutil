#include "align_eval_aux.h"
#include "file_utils.h"

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



std::vector<LineIndex> 
build_index(char const* sam_file,
            char * chunk_buffer,
            size_t max_mem,
            size_t max_line,
            size_t (* samline_pos)(char const*, CONTIG_OFFSETS const&), 
            CONTIG_OFFSETS const& contig_offsets,
            size_t * num_chunks)
{

    std::vector<LineIndex> line_index;
    ssize_t line_length;
    size_t index;

    FILE * sam_fh = fopen(sam_file, "r");
    SetToFirstDataLine(&sam_fh);
    size_t header_length = ftell(sam_fh);
    
    gzFile sam_zfh = gzopen(sam_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(sam_zfh, max_mem, max_line);
    gzclose(sam_zfh);

    chunk_lengths[0] -= header_length;
    *num_chunks = chunk_lengths.size();

    size_t current_offset = header_length;
    for (size_t chunk = 0; chunk != *num_chunks; ++chunk)
    {
        size_t nbytes_read = fread(chunk_buffer, 1, chunk_lengths[chunk], sam_fh);
        assert(nbytes_read == chunk_lengths[chunk]);
        chunk_buffer[nbytes_read] = '\0';

        std::vector<char *> samlines = FileUtils::find_line_starts(chunk_buffer);

        char * line_end = chunk_buffer;
        while ((line_end = strchr(line_end, '\n')) != NULL)
        {
            *line_end++ = '\0';
        }

        std::vector<char *>::iterator sit;
        for (sit = samlines.begin(); sit != samlines.end(); ++sit)
        {
            //zero terminated.  we still want to record the spacer for what was the \n
            line_length = strlen(*sit) + 1; 
            index = samline_pos(*sit, contig_offsets);
            line_index.push_back(LineIndex(index, current_offset, line_length));
            current_offset += line_length;
        }

    }
    fclose(sam_fh);

    return line_index;
}
