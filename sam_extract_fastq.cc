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

#include <omp.h>

#include <parallel/algorithm>

#include "align_eval_raw.h"
#include "align_eval_aux.h"

#include "file_utils.h"
#include "dep/tools.h"
#include "sam_helper.h"
#include "sam_order.h"


int sam_extract_usage(size_t mdef, size_t zdef, size_t ydef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil extract [OPTIONS] input.fsort.sam reads1.fastq.gz [reads2.fastq.gz orphans.fastq.gz]\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-C  STRING    work in the directory named here [.]\n"
            "-z  INT       zlib compression level to be used on tmp files (0-9). 0 means no compression. [%Zu]\n"
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

    char const* tmp_file_prefix = NULL;

    char c;
    while ((c = getopt(argc, argv, "m:t:z:y:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'z': tmp_file_gz_level = static_cast<size_t>(atof(optarg)); break;
        case 'y': tmp_file_gz_strategy = static_cast<size_t>(atof(optarg)); break;
        default: return sam_sort_usage(max_mem_def, tmp_file_gz_level_def, tmp_file_gz_strategy_def); break;
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    __gnu_parallel::_Settings psettings;
    psettings.for_each_minimal_n = 2;
    __gnu_parallel::_Settings::set(psettings);

    size_t argcount = argc - optind;
    bool paired_read_mode = (argcount == 3);

    if (argcount != 1 && argcount != 3)
    {
        return sam_extract_usage(max_mem_def, tmp_file_gz_level_def,
                                 tmp_file_gz_strategy_def);
    }

    char const* sorted_sam_file = argv[optind];
    char const* fastq1_file = argv[optind + 1];

    char const* fastq2_file = "/dev/null";
    char const* orphan_fastq_file = "/dev/null";

    
    if (argcount == 3)
    {
        fastq2_file = argv[optind + 2];
        orphan_fastq_file = argv[optind + 3];
    }
    
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "r");
    FILE * fastq1_fh = open_if_present(fastq1_file, "w");

    // these will be NULL if we are in single read mode
    FILE * fastq2_fh;
    FILE * orphan_fastq_fh;

    if (paired_read_mode)
    {    
        fastq2_fh = open_if_present(fastq2_file, "w");
        orphan_fastq_fh = open_if_present(orphan_fastq_file, "w");
    }

    size_t chunk_size = max_mem / 2;

    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out[2];

    chunk_buffer_out[0] = new char[chunk_size + 1];
    if (paired_read_mode)
    {
        chunk_buffer_out[1] = new char[chunk_size + 1];
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
    SamOrder sam_order(SAM_RID, "FRAGMENT");

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

        clock_gettime(CLOCK_REALTIME, &time_end);
        fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
        fflush(stderr);
            
        if (! sam_order.Initialized() && ! sam_lines.empty())
        {
            sam_order.InitFromID(sam_lines[0]);
        }

        // this won't work if we happen to read exactly the number of
        // bytes left in the file without tripping the 'feof' flag.
        is_last_chunk = feof(sorted_sam_fh);

        size_t num_unprocessed_lines =
            convert_chunk_sam_to_fastq(sam_lines,
                                       & sam_order,
                                       chunk_buffer_in, 
                                       chunk_buffer_out[0],
                                       chunk_buffer_out[1],
                                       orphan_buffer_out,
                                       is_last_chunk);
            

        parallel_compress(chunk_buffer_out[0], chunk_size, num_threads, fastq1_fh);
        if (paired_end_mode)
        {
            parallel_compress(chunk_buffer_out[1], chunk_size, num_threads, fastq2_fh);
            parallel_compress(orphan_buffer_out, chunk_size, num_threads, orphan_fastq_fh);
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

    
    fclose(sorted_sam_fh);
    fclose(fastq1_fh);

    close_if_present(fastq2_fh);
    close_if_present(orphan_fastq_fh);
    
    delete [] chunk_buffer_in;
    delete [] chunk_buffer_out[0];
    if (chunk_buffer_out[1] != NULL)
    {
        delete [] chunk_buffer_out[1];
        delete [] orphan_buffer_out;
    }

    return 0;
}
