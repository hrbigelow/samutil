#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <set>
#include <cassert>

#include "dep/tools.h"

struct coord
{
    char contig[100];
    size_t pos;
    bool is_start;
    coord(char const* _contig, size_t _pos, bool _is_start) :
        pos(_pos), is_start(_is_start) { 
        strcpy(contig, _contig);
    }
};


struct less_coord
{
    bool operator()(coord const& a, coord const& b) const {
        int contig_cmp = strcmp(a.contig, b.contig);
        if (contig_cmp == 0)
        {
            return a.pos < b.pos;
        }
        else
        {
            return contig_cmp < 0;
        }
    }
};



void complete_union_interval(FILE * interval_fh, FILE * ids_fh, 
                             std::multiset<coord, less_coord> * positions,
                             std::set<size_t> * contained_ids, 
                             size_t * union_id)
{
    coord const& start_coord = *(*positions).begin();
    coord const& end_coord = *(*positions).rbegin();

    assert(start_coord.is_start);
    assert(! end_coord.is_start);
    assert(strcmp(start_coord.contig, end_coord.contig) == 0);

    //we need to close out a previous union interval
    fprintf(interval_fh, "%Zu\t%s\t%Zu\t%Zu\n",
            *union_id, start_coord.contig, 
            start_coord.pos, end_coord.pos);

    for (std::set<size_t>::const_iterator 
             contained_ids_iter = (*contained_ids).begin();
         contained_ids_iter != (*contained_ids).end();
         ++contained_ids_iter)
    {
        fprintf(ids_fh, "%Zu\t%Zu\n", *union_id, *contained_ids_iter);
    }
    ++(*union_id);
    (*positions).clear();
    (*contained_ids).clear();
}


int usage()
{
    fprintf(stderr, "Usage:\n\n"
            "union_regions <input_regions> <union_regions> <union2input_idmap>\n\n"
            "<input_regions_file> contains tab-separated fields: id contig start end\n\n"
            "input_regions_file must be sorted as: sort -k 2,2d -k 3,3g -k 4,4g\n\n"
            "all coordinates are boundary-based: (0,10) denotes first 10 bases\n\n"
            "union_regions has format the union intervals in format:\n"
            "id contig start end\n"
            "union2input_idmap: output of union ids to input ids\n");
    exit(1);
}


int main(int argc, char ** argv)
{

    if (argc != 4)
    {
        return usage();
    }

    char * input_regions_file = argv[1];

    char * interval_file = argv[2];
    char * ids_file = argv[3];

    FILE * input_regions_fh = open_or_die(input_regions_file, "r", "Input regions file");

    FILE * interval_fh = fopen(interval_file, "w");
    FILE * ids_fh = fopen(ids_file, "w");

    if (interval_fh == NULL)
    {
        fprintf(stderr, "Error: couldn't open interval file %s for writing\n",
                interval_file);
        exit(1);
    }

    if (ids_fh == NULL)
    {
        fprintf(stderr, "Error: couldn't open ids file %s for writing\n",
                ids_file);
        exit(1);
    }
    

    size_t id;
    char contig[100];
    size_t start, end;

    size_t union_id = 0;
    std::set<size_t> contained_ids;

    std::multiset<coord, less_coord> positions;

    less_coord less_coord_object;

    //assume input is sorted by contig, start, end
    while (fscanf(input_regions_fh, "%zu\t%s\t%zu\t%zu\n", &id, contig, &start, &end) == 4)
    {
        if (end < start)
        {
            fprintf(stderr, "Error: end coordinate is less than start.  Invalid input:\n"
                    "%Zu\t%s\t%Zu\t%Zu\n", id, contig, start, end);
            exit(1);
        }
        coord next_start(contig, start, true);
        coord next_end(contig, end, false);

        if (positions.empty() ||
            less_coord_object(*positions.rbegin(), next_start))
        {
            //we are beginning a new union interval
            if (! positions.empty())
            {
                complete_union_interval(interval_fh, ids_fh, &positions, 
                                        &contained_ids, &union_id);
            }
        }
        assert(positions.empty() || strcmp((*positions.rbegin()).contig, contig) == 0);
        positions.insert(next_start);
        positions.insert(next_end);
        contained_ids.insert(id);
        
    }
    //by definition, there is a gap at the end of the input
    complete_union_interval(interval_fh, ids_fh, &positions, 
                            &contained_ids, &union_id);
    fclose(interval_fh);
    fclose(ids_fh);

}
