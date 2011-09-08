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

#include "align_eval_sort.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <functional>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
//#include <ext/algorithm>


#include "align_eval_raw.h"
#include "align_eval_aux.h"

#include "file_utils.h"
#include "dep/tools.h"
#include "sam_helper.h"
#include "sam_order.h"


std::vector<INDEX_ITER> 
get_quantiles(std::vector<LineIndex> * line_index,
              bool (less_fcn)(LineIndex const&, LineIndex const&),
              size_t num_chunks,
              bool already_sorted)
{
    
    size_t lines_per_chunk = (*line_index).size() / num_chunks;
    std::vector<INDEX_ITER> quantiles;
    
    INDEX_ITER iter = (*line_index).begin();
    INDEX_ITER end = (*line_index).end();

    quantiles.push_back(iter);
    for (size_t n = 0; n != num_chunks - 1; ++n)
    {
        std::advance(iter, lines_per_chunk);
        if (already_sorted)
        {
            //iter already defines the partition.  do nothing.
        }
        else
        {
            std::nth_element(quantiles.back(), iter, end, less_fcn);
        }
        quantiles.push_back(iter);
    }
    quantiles.push_back(end);
    return quantiles;
}



int align_eval_sort_usage(char const* sdef, size_t mdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval sort [OPTIONS] alignment.sam alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-s  STRING    type of sorting to use {READ_ID_FLAG, ALIGN, GUIDE, MIN_ALIGN_GUIDE}[%s]\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-u  FLAG      (unique) if present, omit printing of all but one duplicate lines. [false]\n"
            "-h  STRING    optional sam header if alignment.sam header lacks SQ fields.\n"
            "              If provided, any header lines in alignment.sam will be ignored.\n\n",
            sdef, mdef);

    fprintf(stderr,
            "Sort orders are:\n"
            "READ_ID_FLAG: sort by read id / pair flag (uniquely identifies the physical fragment)\n"
            "ALIGN: sort by alignment position\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_align_eval_sort(int argc, char ** argv)
{

    char const* sort_type_def = "MIN_ALIGN_GUIDE";
    char const* sort_type = sort_type_def;
    char const* sam_header_file = "/dev/null";

    bool filter_duplicates = false;

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    char c;
    while ((c = getopt(argc, argv, "s:m:uh:")) >= 0)
    {
        switch(c)
        {
        case 's': sort_type = optarg; break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'u': filter_duplicates = true; break;
        case 'h': sam_header_file = optarg; break;
        default: return align_eval_sort_usage(sort_type_def, max_mem_def); break;
        }
    }

    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return align_eval_sort_usage(sort_type_def, max_mem_def);
    }

    char const* alignment_sam_file = argv[optind];
    char const* sorted_sam_file = argv[optind + 1];

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
    delete header_buf;

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


    size_t num_chunks;
    char * chunk_buffer = new char[max_mem];
    char * chunk_buffer_out;
    char * write_pointer;

    size_t const max_line = 10000;

    fprintf(stderr, "Building index...");
    fflush(stderr);
    std::vector<LineIndex> line_index = 
        build_index(alignment_sam_file, chunk_buffer, max_mem, max_line,
                    sam_order, &num_chunks);
    fprintf(stderr, "done\n");
    fflush(stderr);

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
    }
    
    // 0. Determine N quantiles based on file offset, offset_quantiles
    // line_index is sorted by (O,O)
    bool already_sorted = true;
    std::vector<INDEX_ITER> offset_quantiles = 
        get_quantiles(&line_index, less_offset, num_chunks, already_sorted);

    std::vector<size_t> offset_quantile_sizes;
    for (size_t o = 0; o != num_chunks; ++o)
    {
        INDEX_ITER beg = offset_quantiles[o];
        INDEX_ITER end = offset_quantiles[o+1];
        
        size_t start_file_offset = (*beg).start_offset;
        size_t end_file_offset = 
            (end == line_index.end())
            ? ((line_index.back()).start_offset + (line_index.back()).line_length)
            : (*end).start_offset;

        offset_quantile_sizes.push_back(end_file_offset - start_file_offset);
    }

    size_t max_offset_quantile =
        *std::max_element(offset_quantile_sizes.begin(), offset_quantile_sizes.end());

    delete chunk_buffer;

    chunk_buffer = new char[max_offset_quantile];
    chunk_buffer_out = new char[max_offset_quantile];

    // 1. For each nth quantile, [n, n+1) of the line index:
    fprintf(stderr, "Sorting index, chunks [0-%zu] ", num_chunks - 1);
    fflush(stderr);

    FILE * out_fh;

    for (size_t o = 0; o != num_chunks; ++o)
    {
        INDEX_ITER beg = offset_quantiles[o];
        INDEX_ITER end = offset_quantiles[o+1];
        
        size_t start_file_offset = (*beg).start_offset;

        write_pointer = chunk_buffer_out;

        //a. Read chunk of file from start.file_offset to end.file_offset
        fseek(alignment_sam_fh, start_file_offset, std::ios::beg);
        fread(chunk_buffer, 1, offset_quantile_sizes[o], alignment_sam_fh);
        
        //b. Sort main index partition [n, n+1) on key value.
        std::sort(beg, end, less_key);

        //c. Print out sorted temp file, or sorted sam file, if there is a single chunk
        if (num_chunks > 1)
        {
            int file_des = mkstemp(tmp_files[o]);
            tmp_fhs[o] = fdopen(file_des, "w");
            out_fh = tmp_fhs[o];
        }
        else
        {
            out_fh = sorted_sam_fh;
        }

        for (INDEX_ITER lit = beg; lit != end; ++lit)
        {
            memcpy(write_pointer,
                   chunk_buffer + (*lit).start_offset - start_file_offset,
                   (*lit).line_length);

            write_pointer += (*lit).line_length;

            // fwrite(chunk_buffer + (*lit).start_offset - start_file_offset, 1,
            //        (*lit).line_length, out_fh);
        }
        assert(static_cast<size_t>(std::distance(chunk_buffer_out, write_pointer)) 
               == offset_quantile_sizes[o]);

        fwrite(chunk_buffer_out, 1, offset_quantile_sizes[o], out_fh);

        if (num_chunks > 1)
        {
            fclose(tmp_fhs[o]);
        }

        fprintf(stderr, "%zu ", o);
        fflush(stderr);
    }
    fprintf(stderr, "done\n");
    fflush(stderr);

    //Now ordered as (O,K) (see step b. above)

    if (num_chunks > 1)
    {

        // 2. Find N quantile key values in index.  (O(N))  Now ordered as (K,K).
        std::vector<INDEX_ITER> key_quantiles = 
            get_quantiles(&line_index, less_key, num_chunks, ! already_sorted);


        // 3. Calculate size of bytes in each key_quantile.
        size_t * key_quantile_sizes = new size_t[num_chunks];
        std::vector<LineIndex> key_quantile_linfo;

        for (size_t k = 0; k != num_chunks; ++k)
        {
            INDEX_ITER kit, kit_end = key_quantiles[k+1];

            LineIndex sentinel = (kit_end == line_index.end())
                ? LineIndex(INT64_MAX, INT64_MAX, 0)
                : *kit_end;

            key_quantile_linfo.push_back(sentinel);

            size_t & kq_size = key_quantile_sizes[k];
            kq_size = 0;
            for (kit = key_quantiles[k]; kit != kit_end; ++kit)
            {
                kq_size += (*kit).line_length;
            }
        }
        size_t kq_size_max = 
            *std::max_element(key_quantile_sizes, key_quantile_sizes + num_chunks);

        // 4. Allocate chunk of memory to hold biggest key_quantile
        delete chunk_buffer;
        delete chunk_buffer_out;

        chunk_buffer = new char[kq_size_max];
        chunk_buffer_out = new char[kq_size_max];

        // 5. Re-partition index by offset_quantiles.  Now ordered as (O,.)
        get_quantiles(&line_index, less_offset, num_chunks, ! already_sorted);

        // order each quantile by key
        for (size_t o = 0; o != num_chunks; ++o)
        {
            std::sort(offset_quantiles[o], offset_quantiles[o+1], less_key);
        }

        // 6. For each key_quantile kq
        for (size_t o = 0; o != num_chunks; ++o)
        {
            tmp_fhs[o] = fopen(tmp_files[o], "r");
        }

        std::vector<INDEX_ITER> prev_sub_k(offset_quantiles.size() - 1);
        std::copy(offset_quantiles.begin(), offset_quantiles.end() - 1, prev_sub_k.begin());

        size_t S, buffer_pos, last_merged_size;
        INDEX_ITER beg, end, sub_k;
        size_t part;
        std::vector<LineIndex> merged;

        size_t num_total_lines = 0;
        size_t num_chunk_lines = 0;

        fprintf(stderr, "Merging chunks [0-%zu] ", num_chunks - 1);
        fflush(stderr);

        for (size_t k = 0; k != num_chunks; ++k)
        {
            write_pointer = chunk_buffer;
            merged.clear();

            //load kq range of each tmp file into the buffer
            buffer_pos = 0;
            for (size_t o = 0; o != num_chunks; ++o)
            {
                beg = prev_sub_k[o];
                end = offset_quantiles[o+1];
                sub_k = std::lower_bound(beg, end, key_quantile_linfo[k], less_key);

                //assert(__gnu_cxx::is_sorted(beg, end, less_key));

                part = std::distance(beg, sub_k);

                // fprintf(stdout, "key %Zu, offset %Zu, num_lines: %Zu\n",
                //         k, o, part);

                S = 0;
                for (INDEX_ITER kit = beg; kit != sub_k; ++kit)
                {
                    (*kit).start_offset = S + buffer_pos;
                    S += (*kit).line_length;
                }
                fread(write_pointer, 1, S, tmp_fhs[o]);
                write_pointer += S;
                buffer_pos += S;

                if (buffer_pos > kq_size_max)
                {
                    fprintf(stderr, "buffer_pos: %Zu, kq_size_max: %Zu, o: %Zu, k: %Zu\n",
                            buffer_pos, kq_size_max, o, k);
                }
                assert(buffer_pos <= kq_size_max);
            

                last_merged_size = merged.size();
                merged.insert(merged.end(), beg, sub_k);

                std::inplace_merge(merged.begin(), merged.begin() + last_merged_size, 
                                   merged.end(), less_key);

                // fprintf(stdout, "key %Zu, offset %Zu, added: %Zu, num_lines: %Zu\n",
                //         k, o, part, merged.size());

                prev_sub_k[o] = sub_k;
            }

            num_total_lines += merged.size();
            num_chunk_lines = merged.size();

            write_pointer = chunk_buffer_out;

            for (INDEX_ITER mit = merged.begin(); mit != merged.end(); ++mit)
            {
                memcpy(write_pointer, chunk_buffer + (*mit).start_offset,
                       (*mit).line_length);

                write_pointer += (*mit).line_length;

                // fwrite(chunk_buffer + (*mit).start_offset, 1,
                //        (*mit).line_length, sorted_sam_fh);
            }

            assert(static_cast<size_t>(std::distance(chunk_buffer_out, write_pointer)) 
                   == key_quantile_sizes[k]);

            fwrite(chunk_buffer_out, 1, key_quantile_sizes[k], sorted_sam_fh);

            fprintf(stderr, "%zu ", k);
            fflush(stderr);

        }

        fprintf(stderr, "done\n");
        fprintf(stderr, "Cleaning up...");
        fflush(stderr);

        delete key_quantile_sizes;
    
        for (size_t o = 0; o != num_chunks; ++o)
        {
            fclose(tmp_fhs[o]);
            remove(tmp_files[o]);
        }
        fprintf(stderr, "done\n");
        fflush(stderr);
    }


    close_if_present(alignment_sam_fh);
    close_if_present(sorted_sam_fh);
    close_if_present(sam_header_fh);

    delete chunk_buffer;
    delete tmp_file_template;
    delete tmp_files;
    delete tmp_fhs;
    delete tmp_file_buf;

    fprintf(stdout, "%Zu lines printed.\n", line_index.size());

    return 0;
}
