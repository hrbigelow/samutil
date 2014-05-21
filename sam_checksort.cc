#include "align_eval_aux.h"
#include "dep/tools.h"
#include "file_utils.h"
#include "sam_index.h"
#include "sam_file.h"

#include <omp.h>

#include <parallel/algorithm>

#include <getopt.h>

int sam_checksort_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "samutil checksort [OPTIONS] sort_order alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-t  INT       number of threads to use [1]\n"
            "-C  STRING    work in the directory named here [.]\n");

    fprintf(stderr,
            "sort_order must be one of:\n"
            "FRAGMENT: sort by read id\n"
            "ALIGN: sort by alignment position\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_sam_checksort(int argc, char ** argv)
{

    size_t num_threads = 1;
    char const* working_dir = ".";

    char c;
    while ((c = getopt(argc, argv, "t:C:")) >= 0)
    {
        switch(c)
        {
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        default: return sam_checksort_usage(); break;
        }
    }

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return sam_checksort_usage();
    }

    char const* sort_type = argv[optind];
    char const* sorted_sam_file = argv[optind + 1];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "r");
    
    contig_dict contig_dictionary;
    index_dict_t flowcell_dict;

    init_contig_length_offset(NULL, & contig_dictionary);

    SAM_INDEX_TYPE index_type = sam_index_type_from_string(sort_type);

    char * header_buf = ReadAllocSAMHeader(sorted_sam_fh);
    delete [] header_buf;

    SetToFirstDataLine(&sorted_sam_fh);

    size_t const chunk_size = 1024l * 1024l * 1024l;

    char * chunk_buffer_in = new char[chunk_size + 1];

    bool is_sorted = true; // innocent until proven guilty
    size_t total_lines = 0;
    size_t nbytes_unused = 0;
    size_t nbytes_read = 0;
    char * last_fragment;
    char * read_pointer = chunk_buffer_in;

    std::vector<char *> samlines;

    while (! feof(sorted_sam_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, sorted_sam_fh);
        
        read_pointer[nbytes_read] = '\0';

        samlines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);
        
        SAM_QNAME_FORMAT qfmt = qname_format(samlines[0]);
        
        sam_index * line_index_chunk = new sam_index[samlines.size()];
        samlines_to_index(num_threads, samlines.data(), samlines.size(), 
                          index_type, qfmt, & contig_dictionary, line_index_chunk, & flowcell_dict);
        
        is_sorted = is_sorted &&
            std::is_sorted(line_index_chunk, line_index_chunk + samlines.size(), samidx_less_key);

        if (! is_sorted)
        {
            for (size_t l = 1; l != samlines.size(); ++l)
            {
                if (samidx_less_key(line_index_chunk[l], line_index_chunk[l-1]))
                // if (line_index_chunk[l].index < line_index_chunk[l-1].index)
                {
                    fprintf(stderr, "Lines %zu and %zu are not in %s order:\n%s\n%s\n\n",
                            total_lines + l - 1,
                            total_lines + l,
                            sort_type,
                            samlines[l - 1],
                            samlines[l]);
                    
                    exit(1);
                }
            }
        }
        delete [] line_index_chunk;
        total_lines += samlines.size();

        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
        
    }
    
    delete [] chunk_buffer_in;

    if (is_sorted)
    {
        fprintf(stderr, "Sorted by %s\n", sort_type);
        return 0;
    }
    else
    {
        fprintf(stderr, "Not sorted by %s\n", sort_type);
        return 1;
    }
}
