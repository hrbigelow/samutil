#include <cstdio>
#include <vector>
#include <parallel/algorithm>

#include "sam_order.h"

int main(int argc, char ** argv)
{
    char * file = argv[1];
    char * nickname = argv[2];
    size_t num_threads = static_cast<size_t>(atof(argv[3]));

    FILE * ids_fh = fopen(file, "r");
    std::vector<size_t> index;
    index.reserve(40000000);
    char qname[1024];

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);

    size_t (* parse_id)(char const* qname);
    char const* format;

    fscanf(ids_fh, "%s\n", qname + 1);
    switch(QNAMEFormat(qname + 1))
    {
    case SAM_NON_INTERPRETED: parse_id = &parse_fragment_id_zero; format = "SAM_NON_INTERPRETED"; break;
    case SAM_NUMERIC: parse_id = &parse_fragment_id_numeric; format = "SAM_NUMERIC"; break;
    case SAM_ILLUMINA: parse_id = &parse_fragment_id_illumina; format = "SAM_ILLUMINA"; break;
    case SAM_CASAVA18: parse_id = &parse_fragment_id_casava_1_8; format = "SAM_CASAVA18"; break;
        
    }
    index.push_back(parse_id(qname + 1));

    while (! feof(ids_fh))
    {
        fscanf(ids_fh, "%s\n", qname);
        index.push_back(parse_id(qname + 1));
    }

    fclose(ids_fh);

    __gnu_parallel::sort(index.begin(), index.end());

    size_t prev_index = 0;
    size_t collisions = 0;

    for (size_t i = 0; i != index.size(); ++i)
    {
        if (! (prev_index < index[i]))
        {
            ++ collisions;
        }
        prev_index = index[i];
    }
    fprintf(stdout, "%s\t%s\t%zu\t%zu\n", nickname, format, index.size(), collisions);
    return 0;
}
