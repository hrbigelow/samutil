//sort an evaluation read by earliest boundary.

/*
  overview:

1. build the index in one file, with file offsets.  keep in memory
2. for each n of N chunks:

   a. partition the index by index value.
   b. sort the nth partition by file offset.
   c. scan through the file, (using litestream) outputting each of the lines
*/


//Maintain a set of 'next' pointers, one for each block.
//Initialize std::vector<char const*> merged to
//size_t m = 0

#include "sam_sort.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <functional>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>

#include <omp.h>

#include <parallel/algorithm>

#include "align_eval_raw.h"
#include "align_eval_aux.h"

#include "file_utils.h"
#include "dep/tools.h"
#include "sam_helper.h"
#include "sam_order.h"


int sam_sort_usage(char const* sdef, size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil sort [OPTIONS] alignment.sam alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-s  STRING    type of sorting to use {READ_ID_FLAG, ALIGN, GUIDE, MIN_ALIGN_GUIDE}[%s]\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-u  FLAG      (unique) if present, omit printing of all but one duplicate lines. [false]\n"
            "-h  STRING    optional sam header if alignment.sam header lacks SQ fields.\n"
            "              If provided, any header lines in alignment.sam will be ignored.\n"
            "-C  STRING    work in the directory named here [.]\n",
            sdef, mdef);

    fprintf(stderr,
            "Sort orders are:\n"
            "READ_ID_FLAG: sort by read id / pair flag (uniquely identifies the physical fragment)\n"
            "ALIGN: sort by alignment position\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_sam_sort(int argc, char ** argv)
{

    char const* sort_type_def = "MIN_ALIGN_GUIDE";
    char const* sort_type = sort_type_def;
    char const* sam_header_file = "/dev/null";

    bool filter_duplicates = false;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;
    char const* working_dir = ".";

    char c;
    while ((c = getopt(argc, argv, "s:m:t:uh:C:")) >= 0)
    {
        switch(c)
        {
        case 's': sort_type = optarg; break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'u': filter_duplicates = true; break;
        case 'h': sam_header_file = optarg; break;
        case 'C': working_dir = optarg; break;
        default: return sam_sort_usage(sort_type_def, max_mem_def); break;
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);


    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return sam_sort_usage(sort_type_def, max_mem_def);
    }

    char const* alignment_sam_file = argv[optind];
    char const* sorted_sam_file = argv[optind + 1];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * alignment_sam_fh = open_if_present(alignment_sam_file, "r");
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "w");
    FILE * sam_header_fh = open_if_present(sam_header_file, "r");

    FILE * used_header_fh = (sam_header_fh != NULL) ? sam_header_fh : alignment_sam_fh;

    SamOrder sam_order(SAM_RID, sort_type);

    sam_order.InitFromFile(alignment_sam_fh);
    sam_order.AddHeaderContigStats(used_header_fh);

    SetToFirstDataLine(&used_header_fh);

    size_t header_length = ftell(used_header_fh);
    rewind(used_header_fh);

    char * header_buf = new char[header_length];
    fread(header_buf, 1, header_length, used_header_fh);
    fwrite(header_buf, 1, header_length, sorted_sam_fh);
    delete [] header_buf;

    fflush(sorted_sam_fh);

    SetToFirstDataLine(&alignment_sam_fh);


    /*
      The index will transition between four possible orderings during these
      steps:
      (O,O) (Partitioned by offset, sub-sorted by offset (i.e. fully sorted by offset)
      (K,K) (Partitioned by key, sub-sorted by key (i.e. fully sorted by key)
      (O,K) (Partitioned by offset, sub-sorted by key)
      (K,O) (Partitioned by key, sub-sorted by offset)

      There are 'num_chunks' (N) partitions.  Each partition is achieved in O(L) time by
      N applications of nth_element(), i.e. O(LN)
      
     */

    /*
      00. Index is initially sorted by global file_offset. (O,O)
      0. Determine N quantiles based on file offset, offset_quantiles
      1. For each nth quantile, [n, n+1) of the line index:
         a. Read chunk of file from start.file_offset to end.file_offset
         b. Sort main index partition [n, n+1) on key value
         c. Print out sorted temp file

      #  Index is now sorted by key in its main partitions.
      2. Find N quantile key values in index.  (O(N))  (now partitioned by key_quantiles)
      3. Calculate size of file bytes in each key_quantile
      4. Allocate chunk of memory to hold biggest key_quantile
      5. Re-partition index by offset_quantiles (O(N)) 

      6. For each key_quantile kq
         aa. Set write_pointer = buffer;
         ab. Create a new, empty merge index called 'merged'
         a. For each nth index partition [n-1, n)
            1. Find sub-partition [q-1, q) in index
            2. Calculate file offset f-1, and size S to load, from sorted index
            3. read S bytes from tmp_file[n] at f-1 into write_pointer. 
            4. write_pointer += S
            5. Update offsets in index at [q-1, q) to be relative to buffer start
               a. index is already sorted by key value in [q-1, q) range
               b. q-1 corresponds to S
               c. each next item is (q-1).line_length after it.
            6. merge(merged.begin(), merged.end(), q-1, q);
            
         b. Print out elements in buffer according to merged, to final sorted file
     */

    size_t const max_line = 10000;
    size_t chunk_size = max_mem / 2;

    gzFile alignment_sam_zfh = gzopen(alignment_sam_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(alignment_sam_zfh, chunk_size, max_line);
    gzclose(alignment_sam_zfh);

    size_t num_chunks = chunk_lengths.size();

    //prepare temporary files
    size_t template_length = strlen(sorted_sam_file) + 8;
    char * tmp_file_template = new char[template_length];
    strcpy(tmp_file_template, sorted_sam_file);
    strcat(tmp_file_template, ".XXXXXX");

    char ** tmp_files = new char *[num_chunks];
    FILE ** tmp_fhs = new FILE *[num_chunks];
    char * tmp_file_buf = new char[num_chunks * (template_length + 1)];
    for (size_t c = 0; c != num_chunks; ++c)
    {
        tmp_files[c] = tmp_file_buf + (c * (template_length + 1));
        strcpy(tmp_files[c], tmp_file_template);
        int fdes = mkstemp(tmp_files[c]);
        tmp_fhs[c] = fdopen(fdes, "w+");
        if (tmp_fhs[c] == NULL)
        {
            fprintf(stderr, "Error: couldn't open temporary chunk file %s for reading/writing.\n",
                    tmp_files[c]);
            exit(1);
        }
    }

    // build the index, and output the
    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out = new char[chunk_size + 1];


    chunk_lengths[0] -= header_length;

    std::vector<size_t> chunk_num_lines;
    std::vector<size_t> offset_quantile_sizes;
    std::vector<LineIndex> line_index;

    fprintf(stderr, "Writing chunks [0-%zu]:", num_chunks - 1);
    fflush(stderr);

    for (size_t c = 0; c != num_chunks; ++c)
    {
        FILE * used_out_fh = num_chunks == 1 ? sorted_sam_fh : tmp_fhs[c];

        size_t nbytes_read = fread(chunk_buffer_in, 1, chunk_lengths[c], alignment_sam_fh);
        assert(nbytes_read == chunk_lengths[c]);
        chunk_buffer_in[nbytes_read] = '\0';

        std::pair<size_t, size_t> chunk_info = 
            process_chunk(chunk_buffer_in, chunk_buffer_out,
                          chunk_lengths[c], sam_order,
                          used_out_fh, & line_index);

        offset_quantile_sizes.push_back(chunk_info.first);
        chunk_num_lines.push_back(chunk_info.second);

        fprintf(stderr, " %zu", c);
        fflush(stderr);
    }

    fprintf(stderr, ".\n");
    fflush(stderr);

    delete [] chunk_buffer_in;
    delete [] chunk_buffer_out;


    // 0. Determine N quantiles based on file offset, offset_quantiles
    std::vector<INDEX_ITER> offset_quantiles;
    INDEX_ITER iit = line_index.begin();
    offset_quantiles.push_back(iit);
    for (size_t c = 0; c != num_chunks; ++c)
    {
        std::advance(iit, chunk_num_lines[c]);
        offset_quantiles.push_back(iit);
    }
    assert(iit == line_index.end());

    // line_index is sorted by (O,K)
    if (num_chunks > 1)
    {
        for (size_t c = 0; c != num_chunks; ++c)
        {
            fseek(tmp_fhs[c], 0, std::ios::beg);
        }

        write_final_merge(line_index, offset_quantiles, tmp_fhs, num_chunks, sorted_sam_fh, NULL);

        fprintf(stderr, ".\n");
        fprintf(stderr, "Cleaning up...");
        fflush(stderr);

        for (size_t o = 0; o != num_chunks; ++o)
        {
            fclose(tmp_fhs[o]);
            remove(tmp_files[o]);
        }
        fprintf(stderr, "done.\n");
        fflush(stderr);
    }


    close_if_present(alignment_sam_fh);
    close_if_present(sorted_sam_fh);
    close_if_present(sam_header_fh);

    delete [] tmp_file_template;
    delete [] tmp_files;
    delete [] tmp_fhs;
    delete [] tmp_file_buf;

    fprintf(stdout, "%Zu lines printed.\n", line_index.size());

    return 0;
}
