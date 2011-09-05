#ifndef _CIGAR_OPS_H
#define _CIGAR_OPS_H

#include <vector>
#include <string>
#include <stdint.h>
#include <set>

#include "seq_projection.h"

//class for applying coordinate transformation implied by CIGAR strings
//

/*

  op   REF   TEMP   SEQ    interpretation

  P      0      0     0    n used only for display purposes

  -      0      0     n    NOT IMPLEMENTED (impossible to generate SEQ from no template)

  H      0      n     0    non-aligned portion of template / read, restricted to ends of CIGAR

  I      0      n     n    insertion to the reference
  S      0      n     n    non-aligned portion of template / read, restricted to ends of CIGAR

  D      n      0     0    deletion from reference due to mutation
  N      n      0     0    deletion from reference due to splicing

  -      n      0     n    NOT IMPLEMENTED (impossible to generate SEQ from no template)

  T      n      n     0    unsequenced (or re-sequenced) bases align
  U      0      n     0    unsequenced (or re-sequenced) bases residing in a section of template that is inserted
                           relative to the reference.

  X      n      n     n    bases align and mismatch
  =      n      n     n    bases align and match
  M      n      n     n    bases align
          
*/


namespace Cigar
{

    enum OpCode
    {
        M = 0, I, D, N, S, H, P, T, U
    };

    //describes which three spans the length operator applies
    struct Op
    {
        OpCode code;
        char name;
        int ref; //whether the op takes up length in the reference
        int temp; 
        int seq;
    };

    extern char const* OpName;

    extern Op Ops[10];

    struct Unit
    {
        Cigar::Op op;
        int64_t length;
        Unit(Cigar::Op const& _op, size_t _l) : op(_op), length(_l) {}
    Unit() : op(Cigar::Ops[Cigar::M]), length(0) { }

    };


    typedef std::vector<Unit> CIGAR_VEC;
    typedef std::vector<Unit>::const_iterator CIGAR_ITER;

    CIGAR_VEC FromString(char const* cigar, int64_t top_to_bottom_offset);

    void ToString(CIGAR_VEC::const_iterator cigar_begin,
                  CIGAR_VEC::const_iterator cigar_end,
                  char * cigar_string);

    std::string ToString(CIGAR_VEC const& cigar);

    size_t UnitLength(Unit const& unit, bool use_top_coord);

    int64_t UnitOffset(Unit const& unit, bool top_to_bottom);

    size_t ReadLength(CIGAR_VEC const& cigar);

    CIGAR_VEC Substring(CIGAR_VEC const& cigar, size_t start_position,
                        size_t end_position, bool use_top_coord);
    
    
    
    size_t Length(CIGAR_VEC const& cigar, bool use_top_coord);
    
    
    size_t Length(CIGAR_ITER cigar_start,
                  CIGAR_ITER cigar_end,
                  bool use_top_coord);

    size_t Overlap(CIGAR_ITER cigar_start, CIGAR_ITER cigar_end);

    //gives the offset from the start of the top sequence to the start
    //of the bottom sequence (if 'top_start_to_bottom_start' is true),
    //otherwise gives the negative of this.
    int64_t LeftOffset(CIGAR_VEC const& cigar, bool top_start_to_bottom_start);


    //gives offset from end of top sequence to end of bottom sequence
    //if 'top_end_to_bottom_end' is true, otherwise the negative of this.
    int64_t RightOffset(CIGAR_VEC const& cigar, bool top_end_to_bottom_end);

    //projects a cigar from (typically) transcriptome space to genome space,
    //using 'blocks' as the projection definition
    CIGAR_VEC Expand(std::vector<block_offsets> const& blocks,
                     CIGAR_VEC const& cigar,
                     bool use_N_as_insert);


    //projects a cigar from (typically) genome space to transcriptome space,
    //using 'blocks' as the projection definition
    CIGAR_VEC Condense(std::vector<block_offsets> const& blocks,
                       CIGAR_VEC const& cigar);

    
    //returns the cigar relationship CIGAR(first, second) from
    //CIGAR(ref, first) and CIGAR(ref, second).
    CIGAR_VEC 
        TransitiveMerge(CIGAR_VEC const& cigar1,
                        std::multiset<std::pair<size_t, size_t> > const& cigar_top_index1,
                        CIGAR_VEC const& cigar2,
                        bool add_padding,
                        bool inserts_are_introns);

    //Used for discovering the agreeing subset of two very similar alignments.
    //inputs: cigar1, (CIGAR(ref, first))
    //        cigar2, (CIGAR(ref,second))
    //output: tcigar: (CIGAR(ref, trimmed_first))
    //'M' states that are the intersection of cigar1 and cigar2
    //are retained in tcigar, the remaining states are converted to 'S'
    CIGAR_VEC 
        TransitiveTrim(CIGAR_VEC const& cigar1,
                       CIGAR_VEC const& cigar2,
                       bool add_padding,
                       bool inserts_are_introns);
        
    
    size_t ProjectCoord(CIGAR_VEC const& transformation,
                        std::multiset<std::pair<size_t, size_t> > const& trans_index,
                        size_t source_coord,
                        bool * projection_applied);

    //project a single coordinate through a CIGAR projection,
    //either from top-to-bottom or bottom-to-top direction
    int64_t ProjectCoord(CIGAR_ITER projection_start,
                         CIGAR_ITER projection_end,
                         int64_t coordinate,
                         bool top_to_bottom,
                         bool * is_missing_projection);
    
    CIGAR_VEC Consolidate(CIGAR_VEC const& source);

    //remove blocks from the ends of source that have zero-length on the top
    //or bottom sequence (top if 'retain_top_extent' is true)
    //after trimming, Length(trimmed, retain_top_extent) will be the same
    void Trim(Cigar::CIGAR_VEC const& source, 
              bool retain_top_extent,
              CIGAR_ITER * trimmed_start,
              CIGAR_ITER * trimmed_end);
    

    CIGAR_VEC Trim(Cigar::CIGAR_VEC const& source, bool retain_top_extent);


    //invert a cigar.  
    CIGAR_VEC Invert(Cigar::CIGAR_VEC const& source, size_t position);

    size_t CountAlignedPositions(Cigar::CIGAR_VEC const& merged,
                                 size_t * total_bases_top,
                                 size_t * total_bases_bottom);

    typedef std::multiset<std::pair<size_t, size_t> > CIGAR_INDEX;
    std::multiset<std::pair<size_t, size_t> > ComputeOffsets(CIGAR_VEC const& cigar);
    
    
}; // end namespace Cigar

#endif // _CIGAR_OPS_H
