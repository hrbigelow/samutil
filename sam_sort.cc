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


int sam_sort_usage(size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil sort [OPTIONS] sort_order alignment.sam alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-u  FLAG      (unique) if present, omit printing of all but one duplicate lines. [false]\n"
            "-h  STRING    optional sam header if alignment.sam header lacks SQ fields.\n"
            "              If provided, any header lines in alignment.sam will be ignored.\n"
            "-g  STRING    optional transcript gtf file.  Required iff sort order is 'PROJALIGN'.\n"
            "-C  STRING    work in the directory named here [.]\n"
            "-T  STRING    path/and/prefix of temp files [name of output file]\n\n",
            mdef);

    fprintf(stderr,
            "sort_order must be one of:\n"
            "FRAGMENT: sort by read id / pair flag (uniquely identifies the physical fragment)\n"
            "ALIGN: sort by alignment position\n"
            "PROJALIGN: sort by genome-projected alignment position.  Requires -g flag\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_sam_sort(int argc, char ** argv)
{

    char const* sam_header_file = "/dev/null";

    bool filter_duplicates = false;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;
    char const* working_dir = ".";

    char const* tmp_file_prefix = NULL;

    char const* gtf_file = NULL;

    char c;
    while ((c = getopt(argc, argv, "m:t:uh:g:C:T:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'u': filter_duplicates = true; break;
        case 'h': sam_header_file = optarg; break;
        case 'g': gtf_file = optarg; break;
        case 'C': working_dir = optarg; break;
        case 'T': tmp_file_prefix = optarg; break;
        default: return sam_sort_usage(max_mem_def); break;
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);


    int arg_count = optind + 3;
    if (argc != arg_count)
    {
        return sam_sort_usage(max_mem_def);
    }

    char const* sort_type = argv[optind];
    char const* alignment_sam_file = argv[optind + 1];
    char const* sorted_sam_file = argv[optind + 2];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    if (tmp_file_prefix == NULL)
    {
        tmp_file_prefix = sorted_sam_file;
    }

    FILE * alignment_sam_fh = open_if_present(alignment_sam_file, "r");
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "w");
    FILE * sam_header_fh = open_if_present(sam_header_file, "r");

    FILE * used_header_fh = (sam_header_fh != NULL) ? sam_header_fh : alignment_sam_fh;

    SamOrder sam_order(SAM_RID, sort_type);
    if (strcmp(sort_type, "PROJALIGN") == 0)
    {
        if (gtf_file == NULL)
        {
            fprintf(stderr, "Error: Sort type given as 'PROJALIGN', but "
                    "no -g option (gtf file) provided\n");
            exit(1);
        }
        else
        {
            sam_order.InitProjection(gtf_file);
        }
    }

    char * header_buf = ReadAllocSAMHeader(used_header_fh);
    size_t header_length = strlen(header_buf);
    fwrite(header_buf, 1, header_length, sorted_sam_fh);

    sam_order.AddHeaderContigStats(header_buf);
    
    delete [] header_buf;

    // merely to slurp up the alignment_sam_fh header
    char * dummy = ReadAllocSAMHeader(alignment_sam_fh);
    delete [] dummy;

    fflush(sorted_sam_fh);

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

    size_t chunk_size = max_mem / 2;

    //prepare temporary file variables
    size_t template_length = strlen(tmp_file_prefix) + 8;
    char * tmp_file_template = new char[template_length];
    strcpy(tmp_file_template, tmp_file_prefix);
    strcat(tmp_file_template, ".XXXXXX");

    std::vector<char *> tmp_files;
    std::vector<FILE *> tmp_fhs;

    // build the index, and output the
    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out = new char[chunk_size + 1];
    char * read_pointer = chunk_buffer_in;

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;

    std::vector<size_t> chunk_num_lines;
    std::vector<size_t> offset_quantile_sizes;
    std::vector<LineIndex> line_index;

    fprintf(stderr, "Writing chunks:");
    fflush(stderr);

    size_t chunk_num = 0;

    while (! feof(alignment_sam_fh))
    {

        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, alignment_sam_fh);
        read_pointer[nbytes_read] = '\0';

        FILE * used_out_fh;
        if (feof(alignment_sam_fh) && tmp_fhs.empty())
        {
            // special case: we can bypass the tmp file stage if the
            // entire input file fits in one chunk
            used_out_fh = sorted_sam_fh;
        }
        else
        {
            char * tmp_file = new char[template_length + 1];
            strcpy(tmp_file, tmp_file_template);
            int fdes = mkstemp(tmp_file);
            FILE * tmp_fh = fdopen(fdes, "w+");
            if (tmp_fh == NULL)
            {
                fprintf(stderr, "Error: couldn't open temporary chunk file %s for reading/writing.\n",
                        tmp_file);
                exit(1);
            }

            tmp_files.push_back(tmp_file);
            tmp_fhs.push_back(tmp_fh);
            
            used_out_fh = tmp_fh;
        }
        
        std::vector<char *> sam_lines = 
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        if (! sam_order.Initialized() && ! sam_lines.empty())
        {
            sam_order.InitFromID(sam_lines[0]);
        }

        std::pair<size_t, size_t> chunk_info = 
            process_chunk(sam_lines, chunk_buffer_in, chunk_buffer_out, 
                          sam_order, used_out_fh, & line_index);

        // now, line_index is in num_chunks pieces, each corresponding
        // to a key-sorted temp file. Each chunk is just a roughly
        // equal portion of the original file and in no particular
        // order. The 'start_offset' fields in line index correspond
        // to original start offsets before key sorting but that are
        // local to the chunk.  (That is, the first entry in each
        // chunk has a start offset of zero.

        offset_quantile_sizes.push_back(chunk_info.first);
        chunk_num_lines.push_back(chunk_info.second);

        fprintf(stderr, " %zu", chunk_num);
        fflush(stderr);

        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
        
        ++chunk_num;
    }

    fprintf(stderr, ".\n");
    fflush(stderr);
    fclose(alignment_sam_fh);

    delete [] chunk_buffer_in;
    delete [] chunk_buffer_out;

    size_t num_chunks = tmp_fhs.empty() ? 1 : tmp_fhs.size();

    // 0. Determine N quantiles based on file offset, offset_quantiles
    set_start_offsets(line_index.begin(), line_index.end(), 0);

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

        write_final_merge(line_index, offset_quantiles, tmp_fhs, sorted_sam_fh, NULL);

        fprintf(stderr, ".\n");
        fprintf(stderr, "Cleaning up...");
        fflush(stderr);

        for (size_t o = 0; o != num_chunks; ++o)
        {
            fclose(tmp_fhs[o]);
            remove(tmp_files[o]);
            delete tmp_files[o];
        }
        fprintf(stderr, "done.\n");
        fflush(stderr);
    }

    close_if_present(sorted_sam_fh);
    close_if_present(sam_header_fh);

    delete [] tmp_file_template;

    fprintf(stderr, "%Zu lines printed.\n", line_index.size());

    return 0;
}
