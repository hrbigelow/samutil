#include "align_eval_aux.h"
#include "dep/tools.h"
#include "sam_order.h"

#include <algorithm>
#include <getopt.h>

int align_eval_checksort_usage(char const* sdef)
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval checksort [OPTIONS] alignment_sorted.sam\n\n"
            "Options:\n\n"
            "-s  STRING    type of sorting to use {READ_ID_FLAG, ALIGN, GUIDE, MIN_ALIGN_GUIDE}[%s]\n",
            sdef);

    fprintf(stderr,
            "Sort orders are:\n"
            "READ_ID_FLAG: sort by read id\n"
            "ALIGN: sort by alignment position\n"
            "GUIDE: sort by read-id encoded guide alignment position\n"
            "MIN_ALIGN_GUIDE: sort by the minimum of ALIGN or GUIDE\n\n");

    return 1;
}


int main_align_eval_checksort(int argc, char ** argv)
{

    char const* sort_type_def = "MIN_ALIGN_GUIDE";
    char const* sort_type = sort_type_def;

    char c;
    while ((c = getopt(argc, argv, "s:")) >= 0)
    {
        switch(c)
        {
        case 's': sort_type = optarg; break;
        default: return align_eval_checksort_usage(sort_type_def); break;
        }
    }

    int arg_count = optind + 1;
    if (argc != arg_count)
    {
        return align_eval_checksort_usage(sort_type_def);
    }

    char const* sorted_sam_file = argv[optind];
    FILE * sorted_sam_fh = open_if_present(sorted_sam_file, "r");
    
    SamOrder sam_order(SAM_RID, sort_type);

    SAM_QNAME_FORMAT qname_fmt = sam_order.InitFromFile(sorted_sam_fh);
    sam_order.AddHeaderContigStats(sorted_sam_fh);

    SamLine::SetGlobalFlags(qname_fmt, "");

    SetToFirstDataLine(&sorted_sam_fh);

    size_t const max_mem = 1024l * 1024l * 100l;
    size_t const max_line = 10000;
    size_t num_chunks;

    char * chunk_buffer = new char[max_mem];

    std::vector<LineIndex> line_index = 
        build_index(sorted_sam_file, chunk_buffer, max_mem, max_line,
                    sam_order, &num_chunks);
    
    delete chunk_buffer;

    bool is_sorted = 
        std::is_sorted(line_index.begin(), line_index.end(), less_key);

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
