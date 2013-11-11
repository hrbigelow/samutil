#include <zlib.h>

#include <parallel/algorithm>

class zstream_tools
{
    size_t tmp_file_gz_strategy;
    size_t tmp_file_gz_level;
    
 public:
    zstream_tools(size_t _strat, size_t _level);

    // compress source_buf using num_threads in parallel, writing compressed output to out_fh
    // returns the number of bytes written
    size_t parallel_compress(char const* source_buf, size_t chunk_size, size_t num_threads, 
                             bool is_last_chunk,
                             FILE * out_fh);
};


struct chunk_buffers {
    unsigned char * source;
    size_t source_len;
    unsigned char * dest;
    size_t dest_len;
    bool last_chunk;
};


// for use with calling zlib's deflate in parallel for_each loop
struct deflate_wrapper
{
    size_t tmp_file_gz_strategy;
    size_t tmp_file_gz_level;
    void operator()(chunk_buffers & cb);
    deflate_wrapper(size_t _strat, size_t _level);
};


