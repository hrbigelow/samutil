#include "align_eval_aux.h"
#include "dep/tools.h"
#include "sam_order.h"
#include "file_utils.h"

#include <omp.h>
#include <zlib.h>

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
    
    SamOrder sam_order(SAM_RID, sort_type);

    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(sorted_sam_fh);
    sam_order.AddHeaderContigStats(sorted_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt, "", "", 0);

    SetToFirstDataLine(&sorted_sam_fh);

    size_t const max_mem = 1024l * 1024l * 1024l;
    size_t const max_line = 4096;

    char * chunk_buffer = new char[max_mem];

    gzFile sorted_sam_zfh = gzopen(sorted_sam_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(sorted_sam_zfh, max_mem, max_line);
    gzclose(sorted_sam_zfh);

    bool is_sorted = true; // innocent until proven guilty
    size_t total_lines = 0;
    size_t nbytes_unused = 0;
    size_t nbytes_read = 0;

    //includes the newline
    size_t last_line_length = 0;

    std::vector<char *> samlines;

    for (size_t c = 0; c != chunk_lengths.size(); ++c)
    {
        nbytes_read = fread(chunk_buffer, 1, chunk_lengths[c], sorted_sam_fh);
        
        chunk_buffer[nbytes_read] = '\0';

        samlines = FileUtils::find_complete_lines_nullify(chunk_buffer, & nbytes_unused);
        assert(nbytes_unused == 0);

        partial_index_aux paux(&sam_order);

        LineIndex * line_index_chunk = new LineIndex[samlines.size()];
        __gnu_parallel::transform(samlines.begin(), samlines.end(), line_index_chunk, paux);
        
        is_sorted = is_sorted &&
            std::is_sorted(line_index_chunk, line_index_chunk + samlines.size(), less_key);

        if (! is_sorted)
        {
            for (size_t l = 1; l != samlines.size(); ++l)
            {
                if (line_index_chunk[l].index < line_index_chunk[l-1].index)
                {
                    char * ln1 = strchr(samlines[l - 1], '\n');
                    *ln1 = '\0';

                    char * ln2 = strchr(samlines[l], '\n');
                    *ln2 = '\0';

                    fprintf(stderr, "Lines %zu and %zu are not in %s order:\n%s\n%s\n\n",
                            total_lines + l - 1,
                            total_lines + l,
                            sort_type,
                            samlines[l - 1],
                            samlines[l]);
                    
                    *ln1 = '\n';
                    *ln2 = '\n';

                    exit(1);
                }
            }
        }
        delete [] line_index_chunk;
        total_lines += samlines.size();
        
        // because we now have to add in a newline
        last_line_length = strlen(samlines.back()) + 1;

        //re-read last line (plus newline)
        if (last_line_length > 0)
        {
            fseek(sorted_sam_fh, -1 * static_cast<off_t>(last_line_length), SEEK_CUR);
        }
    }
    
    delete [] chunk_buffer;

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
