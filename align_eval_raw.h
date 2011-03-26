#ifndef _ALIGN_EVAL_RAW_H
#define _ALIGN_EVAL_RAW_H

#include "align_eval.h"
#include "sam_helper.h"


#include <map>
#include <cstring>

int align_eval_raw_usage(char const* sdef, size_t mdef);
int main_align_eval_raw(int argc, char ** argv);

struct Feature
{
    int guide_jumps : 20;
    int correct_jumps : 20;
    int error_jumps : 24;
    //char : 4;
    OFFSETS_ITER offset_iter;
    Feature() : guide_jumps(0), correct_jumps(0), error_jumps(0) { 
        //memset(this, 0, 8);
    }
    Feature(Feature const& f) : 
        guide_jumps(f.guide_jumps),
        correct_jumps(f.correct_jumps),
        error_jumps(f.error_jumps)
    { }
};



struct read_coords
{
    size_t fragment_id; //unique identifier for the fragment
    size_t read_id; //unique identifier for the read (combination of fragment_id and read-number (1 or 2)
    char contig[1024];
    bool pos_stranded;
    size_t position;
    char cigar[10000];
};


void CigarFromSimSAMLine(char const* readname, bool first_in_pair,
                         bool is_ones_based_position,
                         read_coords * guide_coords);

#endif // _ALIGN_EVAL_RAW_H
