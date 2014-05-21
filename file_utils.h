#ifndef FILE_UTILS_H
#define FILE_UTILS_H


#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <zlib.h>

namespace FileUtils
{
    int switch_printf(bool is_file, void ** buf, char const* format, ...);

    void cat(char * buf, 
             size_t bufsize, 
             size_t bytes_to_print, 
             FILE * input_fh, 
             FILE * output_fh);

    void print_chunk_by_line(FILE * input_fh, size_t chunk, 
                             size_t chunk_approx_size,
                             FILE * output_fh);

    std::vector<size_t> ChunkLengths(gzFile input_fh, size_t chunk_approx_size,
                                     size_t max_line);

    std::vector<size_t> chunk_lengths(FILE * input_fh, size_t chunk_approx_size,
                                      size_t max_line);
    
    std::vector<char *> find_complete_lines(char * lines, char ** last_fragment);

    std::vector<char *> find_complete_lines_nullify(char * lines, char ** last_fragment);

    size_t read_until_newline(char * buffer, size_t bytes_wanted, size_t max_line_length, FILE * fh,
                              size_t * fread_elapsed_nsec);

    void print_until_delim(char const* data, char delim, FILE * out_fh);

};


class BufferedFile
{
    char * chunk_buffer;
    size_t chunk_size;
    size_t max_line;
    gzFile fh;

 public:

    char * file;
    bool valid;

    std::vector<size_t> chunk_lengths;
    std::vector<size_t>::const_iterator chunk_iter;

    std::vector<char *> line_starts;
    std::vector<char *>::iterator line_iter;

    BufferedFile(char const* _file, size_t _chunk_size, size_t _max_line);
    ~BufferedFile();
    bool initialize();

    char * next_n_lines(size_t num_lines, bool * advanced_chunk, bool * reloaded_file);
};

#endif // _FILE_UTILS_H
