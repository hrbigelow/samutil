/*
  In-memory parallel compression using zlib.h
 */

#include "zstream_tools.h"
#include "gzip_tools.h"
#include "time_tools.h"

#include <cassert>

deflate_wrapper::deflate_wrapper(size_t _strat, size_t _level) :
    tmp_file_gz_strategy(_strat),
    tmp_file_gz_level(_level) { }


// ensures that we only process UINT_MAX amounts at a time
void deflate_wrapper::operator()(chunk_buffers & cb)
{
    if (cb.source_len == 0)
    {
        cb.dest_len = 0;
        return;
    }
    if (this->tmp_file_gz_level == 0)
    {
        cb.dest = cb.source;
        cb.dest_len = cb.source_len;
        return;
    }
    
    z_stream zs;
    zs.zalloc = Z_NULL;
    zs.zfree = Z_NULL;
    zs.opaque = Z_NULL;

    // write a raw deflate stream, use max memory for performance
    deflateInit2(&zs, this->tmp_file_gz_level, Z_DEFLATED, -15, 9, this->tmp_file_gz_strategy);

    zs.next_in = cb.source;
    zs.avail_in = 0;
    zs.next_out = cb.dest;
    zs.avail_out = 0;

    // will be compressing a total of 'left' bytes, in chunks of no more than MAXP2 bytes each.
    size_t left = cb.source_len;

    // deflate chunks of MAXP2 bytes at a time.  Use Z_FINISH at the last chunk.
    unsigned int chunk_len = 0;
    unsigned int prev_chunk_len = 0;
    int flush;
    int deflate_retval;

    while (left != 0)
    {
        // update the next_in pointer
        chunk_len = left > MAXP2 ? MAXP2 : left;
        zs.next_in += prev_chunk_len - zs.avail_in;
        zs.avail_in = chunk_len;
        flush = (left > MAXP2) ? Z_NO_FLUSH : (cb.last_chunk ? Z_FINISH : Z_SYNC_FLUSH);
        zs.next_out += prev_chunk_len - zs.avail_out;
        zs.avail_out = chunk_len;
        deflate_retval = deflate(& zs, flush);
        left -= chunk_len;
        prev_chunk_len = chunk_len;
    }
    // assert(zs.avail_in == 0);

    // assert(deflate_retval == Z_STREAM_END);

    // inform the caller how many bytes were written
    cb.dest_len = zs.total_out;

    deflateEnd(&zs);

}


// void deflate_wrapper::operator()(z_stream & zs)
// {
//     // use Z_FINISH since we are only using these zstreams once.
//     int ret = deflate(& zs, Z_FINISH);
//     assert(ret == Z_STREAM_END);
// }


zstream_tools::zstream_tools(size_t _strat, size_t _level) : 
    tmp_file_gz_strategy(_strat), tmp_file_gz_level(_level) { }




// return the number of bytes writte
// splits up the source_buf into num_threads chunks
// 
size_t zstream_tools::parallel_compress(char const* source_buf, 
                                        size_t chunk_size, 
                                        size_t num_threads,
                                        bool is_last_chunk,
                                        FILE * out_fh)
{

    // deflate doesn't work properly if it has zero size, or does it?

    timespec time_begin, time_end;

    std::vector<chunk_buffers> parts(num_threads);

    char * z_chunk_buffer_out = new char[chunk_size];

    size_t subchunk_size = chunk_size / num_threads;
    size_t extra_at_end = chunk_size % num_threads;

    for (size_t t = 0; t != num_threads; ++t)
    {
        parts[t].source = reinterpret_cast<unsigned char *>(const_cast<char *>(source_buf) + (t * subchunk_size));
        parts[t].source_len = (t == num_threads - 1) ? (subchunk_size + extra_at_end) : subchunk_size;
        parts[t].dest = reinterpret_cast<unsigned char *>(z_chunk_buffer_out + (t * subchunk_size));
        parts[t].last_chunk = (is_last_chunk && t == (num_threads - 1));
    }


    clock_gettime(CLOCK_REALTIME, &time_begin);
    fprintf(stderr, "Compressing...");
    fflush(stderr);

    deflate_wrapper dwrap(this->tmp_file_gz_strategy,
                          this->tmp_file_gz_level);

    // parallel deflate, writing the deflated contents into z_chunk_buffer_out
    __gnu_parallel::for_each(parts.begin(), parts.end(), dwrap);

    clock_gettime(CLOCK_REALTIME, &time_end);

    fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
    fflush(stderr);

    clock_gettime(CLOCK_REALTIME, &time_begin);
    fprintf(stderr, "Writing to file...");
    fflush(stderr);

    // sequential write of the subchunks into the tmp_fh
    size_t nbytes_written = 0;
    for (size_t t = 0; t != num_threads; ++t)
    {
        nbytes_written += fwrite(parts[t].dest, 1, parts[t].dest_len, out_fh);
    }

    clock_gettime(CLOCK_REALTIME, &time_end);
    fprintf(stderr, "%Zu bytes written. %Zu ms\n", nbytes_written, elapsed_ms(time_begin, time_end));
    fflush(stderr);

    // clean up
    delete z_chunk_buffer_out;

    return nbytes_written;
}


