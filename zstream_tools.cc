/*
  In-memory parallel compression using zlib.h
 */

#include "zstream_tools.h"
#include "time_tools.h"



void deflate_wrapper::operator()(z_stream & zs)
{
    (void) deflate(& zs, Z_FINISH);
}

zstream_tools::zstream_tools(size_t _strat, size_t _level) : 
    tmp_file_gz_strategy(_strat), tmp_file_gz_level(_level) { }



// return the number of bytes written
size_t zstream_tools::parallel_compress(char const* source_buf, 
                                        size_t chunk_size, 
                                        size_t num_threads, FILE * out_fh)
{

    timespec time_begin, time_end;

    std::vector<z_stream> zstreams(num_threads);
    std::vector<unsigned char *> z_chunk_starts(num_threads);
    std::vector<unsigned int> input_chunk_lengths(num_threads);

    char * z_chunk_buffer_out = new char[chunk_size];

    size_t subchunk_size = chunk_size / num_threads;
    size_t extra_at_end = chunk_size % num_threads;
    for (size_t t = 0; t != num_threads; ++t)
    {
        input_chunk_lengths[t] = (t == num_threads - 1) ? (subchunk_size + extra_at_end) : subchunk_size; 
        z_chunk_starts[t] = (unsigned char *)(z_chunk_buffer_out + (t * subchunk_size));
        zstreams[t].zalloc = Z_NULL;
        zstreams[t].zfree = Z_NULL;
        zstreams[t].opaque = Z_NULL;
        deflateInit2(&zstreams[t], this->tmp_file_gz_level, 8, 15 + 16, 8, this->tmp_file_gz_strategy); // from minigzip.c, for writing gzip
        zstreams[t].next_in = (unsigned char *)(source_buf + (t * subchunk_size));
        zstreams[t].avail_in = input_chunk_lengths[t];
        zstreams[t].next_out = z_chunk_starts[t];
        zstreams[t].avail_out = input_chunk_lengths[t];
    }

    clock_gettime(CLOCK_REALTIME, &time_begin);
    fprintf(stderr, "Compressing...");
    fflush(stderr);

    // parallel deflate, writing the deflated contents into z_chunk_buffer_out
    __gnu_parallel::for_each(zstreams.begin(), zstreams.end(), deflate_wrapper());

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
        nbytes_written += fwrite(z_chunk_starts[t], 1, input_chunk_lengths[t] - zstreams[t].avail_out, out_fh);
        deflateEnd(& zstreams[t]);
    }

    clock_gettime(CLOCK_REALTIME, &time_end);
    fprintf(stderr, "%Zu bytes written. %Zu ms\n", nbytes_written, elapsed_ms(time_begin, time_end));
    fflush(stderr);

    // clean up
    delete z_chunk_buffer_out;

    return nbytes_written;
}
