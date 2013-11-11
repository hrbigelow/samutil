/*
Extract fastq-formatted data from a SAM file.
SAM file must be given in FRAGMENT sorted order.

Reverse-complements the SEQ and QUAL fields based on the 0x10 SAM flag.
Outputs the first seen read of each fragment / read pair combination to appropriate files.

Allows user to specify whether to expect single or paired end reads.

If expecting paired-end reads, outputs all pairs to fastq1 and fastq2
files, and outputs orphans to an 'orphan' file.

Reads are determined to be orphans if NOT both reads in the pair are
encountered.

*/

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <functional>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <zlib.h>
#include <ctime>
#include <time.h>

#include <omp.h>

#include <parallel/algorithm>

#include "sam_to_fastq.h"
#include "time_tools.h"
#include "zstream_tools.h"
#include "gzip_tools.h"

#include "file_utils.h"


int sam_extract_usage(size_t mdef, size_t zdef, size_t ydef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil extract [OPTIONS] input.fsort.sam reads1.fastq.gz [reads2.fastq.gz orphans.fastq.gz]\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-z  INT       zlib compression level to be used on tmp files (0-8). 0 means no compression. [%Zu]\n"
            "-y  INT       zlib compression strategy for tmp files (0-4).  [%Zu]\n"
            "              0=Z_DEFAULT_STRATEGY, 1=Z_FILTERED, 2=Z_HUFFMAN_ONLY, 3=Z_RLE, 4=Z_FIXED\n\n",
            mdef, zdef, ydef);

    fprintf(stderr,
            "If expecting single-end reads, call with reads1.fastq.gz\n"
            "If expecting paired-end reads, call with reads1.fastq.gz reads2.fastq.gz orphans.fastq.gz\n\n");

    return 1;
}


int main_sam_extract(int argc, char ** argv)
{

    size_t tmp_file_gz_level_def = 2;
    size_t tmp_file_gz_level = tmp_file_gz_level_def;

    size_t tmp_file_gz_strategy_def = 0;
    size_t tmp_file_gz_strategy = tmp_file_gz_strategy_def;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;

    char c;
    while ((c = getopt(argc, argv, "m:t:z:y:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'z': tmp_file_gz_level = static_cast<size_t>(atof(optarg)); break;
        case 'y': tmp_file_gz_strategy = static_cast<size_t>(atof(optarg)); break;
        default: return sam_extract_usage(max_mem_def, tmp_file_gz_level_def, tmp_file_gz_strategy_def); break;
        }
    }

    if (tmp_file_gz_level > 8)
    {
        fprintf(stderr, "Error: zlib compression level > 8 is not supported\n");
        return sam_extract_usage(max_mem_def, tmp_file_gz_level_def,
                                 tmp_file_gz_strategy_def);
    }
    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    __gnu_parallel::_Settings psettings;
    psettings.for_each_minimal_n = 2;
    __gnu_parallel::_Settings::set(psettings);

    size_t argcount = argc - optind;
    bool paired_read_mode = (argcount == 4);

    if (argcount != 2 && argcount != 4)
    {
        return sam_extract_usage(max_mem_def, tmp_file_gz_level_def,
                                 tmp_file_gz_strategy_def);
    }

    char const* sorted_sam_file = argv[optind];
    char const* fastq1_file = argv[optind + 1];

    char const* fastq2_file = "/dev/null";
    char const* orphan_fastq_file = "/dev/null";

    
    if (paired_read_mode)
    {
        fastq2_file = argv[optind + 2];
        orphan_fastq_file = argv[optind + 3];
    }
    
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "r");
    FILE * fastq1_fh = open_if_present(fastq1_file, "w");

    // these will be NULL if we are in single read mode
    FILE * fastq2_fh = NULL;
    FILE * orphan_fastq_fh = NULL;


    size_t chunk_size;
    if (paired_read_mode)
    {    
        fastq2_fh = open_if_present(fastq2_file, "w");
        orphan_fastq_fh = open_if_present(orphan_fastq_file, "w");
        chunk_size = max_mem / 4;
    }
    else
    {
        chunk_size = max_mem / 2;
    }

    FILE * out_fhs[3] = { fastq1_fh, fastq2_fh, orphan_fastq_fh };

    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out[2] = { NULL, NULL };
    char * orphan_buffer_out = NULL;

    chunk_buffer_out[0] = new char[chunk_size + 1];
    if (paired_read_mode)
    {
        chunk_buffer_out[1] = new char[chunk_size + 1];
        orphan_buffer_out = new char[chunk_size + 1];
    }

    char * read_pointer = chunk_buffer_in;

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;

    std::vector<size_t> chunk_num_lines;

    size_t chunk_num = 0;
    timespec time_begin, time_end;

    // advance sorted_sam_fh to the reads section.
    char * dummy = ReadAllocSAMHeader(sorted_sam_fh);
    delete dummy;

    bool is_last_chunk;

    zstream_tools zt(tmp_file_gz_strategy, tmp_file_gz_level);
    
    time_t current_time = time(NULL);
    // unsigned long fq_headlen[3];

    if (tmp_file_gz_level > 0)
    {
        put_gzip_header(NULL, current_time, tmp_file_gz_level, fastq1_fh);
    }

    if (paired_read_mode && tmp_file_gz_level > 0)
    {
        put_gzip_header(NULL, current_time, tmp_file_gz_level, fastq2_fh);
        put_gzip_header(NULL, current_time, tmp_file_gz_level, orphan_fastq_fh);
    }

    size_t total_raw_bytes[3] = { 0, 0, 0 };
    unsigned long crc32_check[3] = { crc32(0L, Z_NULL, 0), crc32(0L, Z_NULL, 0), crc32(0L, Z_NULL, 0) };

    while (! feof(sorted_sam_fh))
    {

        clock_gettime(CLOCK_REALTIME, &time_begin);
        fprintf(stderr, "Reading chunk %Zu...", chunk_num + 1);
        fflush(stderr);

        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, sorted_sam_fh);
        clock_gettime(CLOCK_REALTIME, &time_end);
        fprintf(stderr, "%Zu bytes read. %Zu ms\n", nbytes_read, elapsed_ms(time_begin, time_end));
        fflush(stderr);

        read_pointer[nbytes_read] = '\0';

        clock_gettime(CLOCK_REALTIME, &time_begin);
        fprintf(stderr, "nullifying lines...");
        fflush(stderr);

        std::vector<char *> sam_lines = 
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        size_t S = sam_lines.size();
        SAMFastqView * fastq_view = new SAMFastqView[S];

        clock_gettime(CLOCK_REALTIME, &time_end);
        fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
        fflush(stderr);
            
        // clunky test for whether there is anything left to read
        char test_char = getc(sorted_sam_fh);
        ungetc(test_char, sorted_sam_fh);
        is_last_chunk = feof(sorted_sam_fh);

        init_fastq_view(sam_lines, fastq_view);

        size_t s_last = find_last_fragment_bound(fastq_view, S, is_last_chunk);
        size_t num_unprocessed_lines = S - s_last;
        // printf("s_last = %Zu\n", s_last);

        size_t nbytes_written[3];

        if (paired_read_mode)
        {
            clock_gettime(CLOCK_REALTIME, &time_begin);
            fprintf(stderr, "Converting from sam to fastq (in memory)...");
            fflush(stderr);

            set_flags_fastq_view_paired(sam_lines, fastq_view, is_last_chunk);
            write_fastq_view_paired(fastq_view, s_last, chunk_buffer_out[0],
                                    chunk_buffer_out[1], orphan_buffer_out,
                                    nbytes_written);

            clock_gettime(CLOCK_REALTIME, &time_end);
            fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
            fflush(stderr);

            char * bufs[3] = { chunk_buffer_out[0], chunk_buffer_out[1], orphan_buffer_out };

            clock_gettime(CLOCK_REALTIME, &time_begin);
            fprintf(stderr, "Computing CRC32 of fastq chunk...");
            fflush(stderr);

            for (size_t i = 0; i != 3; ++i)
            {
                crc32_check[i] = 
                    crc32_comb(crc32_check[i], 
                               crc32_chunk(bufs[i], nbytes_written[i]),
                               nbytes_written[i]);
                total_raw_bytes[i] += nbytes_written[i];
                // printf("total_raw_bytes[%Zu] += %Zu\n", i, nbytes_written[i]);
            }

            clock_gettime(CLOCK_REALTIME, &time_end);
            fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
            fflush(stderr);

        }
        else
        {
            set_flags_fastq_view_single(sam_lines, fastq_view, is_last_chunk);
            write_fastq_view_single(fastq_view, s_last, chunk_buffer_out[0], nbytes_written);

            crc32_check[0] = 
                crc32_comb(crc32_check[0], 
                           crc32_chunk(chunk_buffer_out[0], nbytes_written[0]),
                           nbytes_written[0]);

            total_raw_bytes[0] += *nbytes_written;
            
        }

        delete [] fastq_view;

        zt.parallel_compress(chunk_buffer_out[0], nbytes_written[0], num_threads, is_last_chunk, fastq1_fh);
        if (paired_read_mode)
        {
            zt.parallel_compress(chunk_buffer_out[1], nbytes_written[1], num_threads, is_last_chunk, fastq2_fh);
            zt.parallel_compress(orphan_buffer_out, nbytes_written[2], num_threads, is_last_chunk, orphan_fastq_fh);
        }

        // now need to put back the unprocessed lines plus number of
        // bytes unused from the last chunk. copy unused samlines into
        // chunk_buffer_in need to de-nullify the lines as well.
        read_pointer = chunk_buffer_in;
        for (size_t s = S - num_unprocessed_lines; s != S; ++s)
        {
            strcpy(read_pointer, sam_lines[s]);
            read_pointer += strlen(sam_lines[s]);
            *read_pointer = '\n';
            ++read_pointer;
        }
        nbytes_unused = strlen(last_fragment);
        memmove(read_pointer, last_fragment, nbytes_unused);
        read_pointer += nbytes_unused;

        nbytes_unused = read_pointer - chunk_buffer_in;
        ++chunk_num;
    }        

    for (size_t f = 0; f != 3; ++f)
    {
        if (out_fhs[f] != NULL && tmp_file_gz_level > 0)
        {
            put_gzip_trailer(total_raw_bytes[f], crc32_check[f], out_fhs[f]);
        }
        if (total_raw_bytes[f] == 0)
        {
            fflush(out_fhs[f]);
            ftruncate(fileno(out_fhs[f]), 0);
        }
        close_if_present(out_fhs[f]);
    }
    
    fclose(sorted_sam_fh);
    
    delete chunk_buffer_in;
    delete chunk_buffer_out[0];
    if (chunk_buffer_out[1] != NULL)
    {
        delete chunk_buffer_out[1];
        delete orphan_buffer_out;
    }

    return 0;
}
