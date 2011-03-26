#ifndef _CIGAR_OPS_H
#define _CIGAR_OPS_H

#include <vector>
#include <string>
#include <stdint.h>
#include <set>

//class for applying coordinate transformation implied by CIGAR strings
//

namespace Cigar
{

    enum Op 
    {
        M = 0, I, D, N, S, H, P, None
        /*
          M: Alignment match (can be a sequence match or mismatch)
          I: Insertion to the reference
          D: Deletion from the reference
          N: Skipped region from the reference
          S: Soft clip on the read (clipped sequence present in <seq>)
          H: Hard clip on the read (clipped sequence NOT present in <seq>)
          P: Padding (silent deletion from the padded reference sequence)

          MP   blocks are positive length for both sequences (reference and read)
          DN   blocks all are positive length for the top sequence (reference)
          I    blocks positive length for bottom sequence (read)
          SH   blocks are zero length for both sequences

          The ops and inverses are:
          ops:       MIDNSHP
          inverses:  MDII??P
         */
    };


    enum OpComp 
    {
        MM = 8, MI, MD, MN, MS, MH, MP,
        IM, II, ID, IN, IS, IH, IP,
        DM, DI, DD, DN, DS, DH, DP,
        NM, NI, ND, NN, NS, NH, NP,
        SM, SI, SD, SN, SS, SH, SP,
        HM, HI, HD, HN, HS, HH, HP,
        PM, PI, PD, PN, PS, PH, PP
    };

    struct Unit
    {
        Cigar::Op op;
        size_t length;
        Unit(Cigar::Op const& _op, size_t _l) : op(_op), length(_l) {}
        Unit() : op(Cigar::M), length(0) { }
        /* Unit(Unit const& u) : op(u.op), length(u.length) { */
        /*     fprintf(stdout, "Unit(Unit const& u)\n"); */
        /* } */
        /* Unit(Unit & u) : op(u.op), length(u.length) {  */
        /*     fprintf(stdout, "Unit(Unit & u)\n"); */
        /* } */
        /* ~Unit() {  */
        /*     fprintf(stdout, "~Unit()\n"); */
        /* } */

        /* Unit & operator=(Unit const& u) { */
        /*     fprintf(stdout, "Unit & operator=(Unit const& u)\n"); */
        /*     this->op = u.op; */
        /*     this->length = u.length; */
        /*     return *this; */
        /* } */


    };


    struct CompareSingle
    {
        Cigar::Op guide_op;
        Cigar::Op test_op;
        int64_t offset;
    };


    struct UnitComparison
    {
        Cigar::OpComp op_pair;
        Cigar::Op guide_op;
        Cigar::Op test_op;
        size_t length;
        int64_t offset;
    };


    Cigar::Op Char2Op(char op);

    char Op2Char(Cigar::Op op);

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

    //returns the cigar relationship CIGAR(first, second) from
    //CIGAR(ref, first) and CIGAR(ref, second).
    CIGAR_VEC 
        TransitiveMerge(CIGAR_VEC const& cigar1,
                        std::multiset<std::pair<size_t, size_t> > const& cigar_top_index1,
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

    size_t CountAlignedPositions(Cigar::CIGAR_VEC const& merged);

    typedef std::multiset<std::pair<size_t, size_t> > CIGAR_INDEX;
    std::multiset<std::pair<size_t, size_t> > ComputeOffsets(CIGAR_VEC const& cigar);
    
    
}; // end namespace Cigar

#endif // _CIGAR_OPS_H
