#include "sam_filter_aux.h"
#include "file_utils.h"
#include "sam_line.h"
#include "sam_file.h"

#include <string.h>
#include <vector>
#include <parallel/algorithm>

extern "C" {
#include "tools.h"
}

int filter_usage(size_t mdef)
{
    fprintf(stderr,
            "\n\nUsage:\n\n"
            "samutil filter [OPTIONS] input.sam filtered.sam\n\n"
            "Options:\n\n"
            "-m  INT       number bytes of memory to use [%Zu]\n"
            "-t  INT       number of threads to use [1]\n"
            "-f  STRING    file containing list of fragment IDs for use with -i []\n"
            "-i  FLAG      (include) if present, only output records with fragment IDs\n"
            "              listed in -f.  Otherwise, only output records with fragment IDs\n"
            "              NOT listed in -f. [absent]\n",
            mdef);

    return 1;
}


int main_sam_filter(int argc, char **argv)
{
    
    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;
    size_t num_threads = 1;

    bool do_include_listed_ids = false;
    char const* fragment_id_file = "/dev/null";

    char c;
    while ((c = getopt(argc, argv, "m:t:f:i")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'f': fragment_id_file = optarg; break;
        case 'i': do_include_listed_ids = true; break;
        default: return filter_usage(max_mem_def); break;
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);
    
    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return filter_usage(max_mem_def);
    }

    char const* in_sam_file = argv[optind];
    char const* out_sam_file = argv[optind + 1];

    FILE * fragment_id_fh = open_if_present(fragment_id_file, "r");
    FILE * in_sam_fh = open_if_present(in_sam_file, "r");
    FILE * out_sam_fh = open_if_present(out_sam_file, "w");

    size_t nbytes_read, nbytes_unused = 0;

    size_t chunk_size = max_mem / 2;

    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out = new char[chunk_size + 1];
    char * read_pointer = chunk_buffer_in;

    // parse the fragment_id_fh
    std::set<std::string> fragment_ids;
    char * id = new char[256];
    ssize_t nfields_read;
    while (! feof(fragment_id_fh))
    {
        nfields_read = fscanf(fragment_id_fh, "%s\n", id);
        fragment_ids.insert(std::string(id));
    }
    fclose(fragment_id_fh);
    delete id;

    bool * ids_found;
    samline_fragment_is_member fragment_ids_check(fragment_ids);
    char * last_fragment;

    char * header = ReadAllocSAMHeader(in_sam_fh);
    size_t header_length = strlen(header);
    SetToFirstDataLine(& in_sam_fh);
    fwrite(header, 1, header_length, out_sam_fh);
    delete[] header;

    while (! feof(in_sam_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, in_sam_fh);
        read_pointer[nbytes_read] = '\0';

        std::vector<char *> sam_lines = 
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        ids_found = new bool[sam_lines.size()];
        std::fill(ids_found, ids_found + sam_lines.size(), false);

        // process the null-terminated SAM lines
        __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), ids_found, fragment_ids_check);

        // iterate through
        size_t line_num = 0;
        char * write_pointer = chunk_buffer_out;
        size_t line_len;
        std::vector<char *>::iterator lit;
        for (lit = sam_lines.begin(); lit != sam_lines.end(); ++lit, ++line_num)
        {
            if (ids_found[line_num] == do_include_listed_ids)
            {
                line_len = strlen(*lit);
                memcpy(write_pointer, *lit, line_len);
                write_pointer[line_len] = '\n';
                write_pointer += line_len + 1;
            }
            else
            {
                // fprintf(stderr, "%s\n", *lit);
            }
        }

        size_t nbytes_written = std::distance(chunk_buffer_out, write_pointer);

        fwrite(chunk_buffer_out, 1, nbytes_written, out_sam_fh);
        fflush(out_sam_fh);

        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

        delete[] ids_found;
    }
    fclose(in_sam_fh);
    fclose(out_sam_fh);

    delete chunk_buffer_in;
    delete chunk_buffer_out;

    return 0;
}
