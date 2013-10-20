#include <zlib.h>

#include <parallel/algorithm>

// for use with calling zlib's deflate in parallel for_each loop
struct deflate_wrapper
{
    void operator()(z_stream & zs);
};

// compress source_buf using num_threads in parallel, writing compressed output to out_fh
void parallel_compress(char const* source_buf, size_t chunk_size, size_t num_threads, FILE * out_fh);

