#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "cigar_ops.h"

#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>


namespace Cigar
{
    char const* OpName = "MIDNSHPTU";

    Op Ops[] = 
        { 
            { M, 'M', 1, 1, 1 },
            { I, 'I', 0, 1, 1 },
            { D, 'D', 1, 0, 0 },
            { N, 'N', 1, 0, 0 },
            { S, 'S', 0, 1, 1 },
            { H, 'H', 0, 1, 0 },
            { P, 'P', 0, 0, 0 },
            { T, 'T', 1, 1, 0 },
            { U, 'U', 0, 1, 0 }
        };
    
};


//parse a cigar string into a vector of cigar units
Cigar::CIGAR_VEC Cigar::FromString(char const* cigar, 
                                   int64_t top_to_bottom_offset)
{
    Cigar::CIGAR_VEC cigar_vec;
    char opname;
    int64_t oplength;

    if (top_to_bottom_offset != 0)
    {
        Cigar::Op first_op = top_to_bottom_offset > 0 ? Cigar::Ops[Cigar::D] : Cigar::Ops[Cigar::I];
        int64_t first_length = labs(top_to_bottom_offset);
        Cigar::Unit unit(first_op, first_length);
        cigar_vec.push_back(unit);
    }
    
    char const* chunk = cigar;
    while (*chunk != '\0')
    {
        size_t num_fields = sscanf(chunk, "%zi%c", &oplength, &opname);
        assert(num_fields == 2);

        size_t op_ind = std::distance(Cigar::OpName, strchr(Cigar::OpName, opname));
        Cigar::Unit unit(Cigar::Ops[op_ind], oplength);

        cigar_vec.push_back(unit);

        chunk = strchr(chunk, opname) + 1;
        assert(chunk != NULL);
    }

    assert(cigar_vec.size() != 0);
    return Cigar::Consolidate(cigar_vec);
}


void Cigar::ToString(Cigar::CIGAR_VEC::const_iterator cigar_start,
                     Cigar::CIGAR_VEC::const_iterator cigar_end,
                     char * cigar_string)
{
    Cigar::CIGAR_VEC::const_iterator iter;
    cigar_string[0] = '\0';
    char op_string[32];
    op_string[1] = '\0';
    
    for (iter = cigar_start; iter != cigar_end; ++iter)
    {
        sprintf(op_string, "%zi%c", (*iter).length, (*iter).op.name);
        strcat(cigar_string, op_string);
    }
}


std::string Cigar::ToString(Cigar::CIGAR_VEC const& cigar)
{
    std::string cigar_string;
    char op_string[100];
    
    for (Cigar::CIGAR_VEC::const_iterator iter = cigar.begin(); 
         iter != cigar.end(); ++iter)
    {
        sprintf(op_string, "%zi%c", (*iter).length, (*iter).op.name);
        cigar_string += op_string;
    }
    return cigar_string;
}


size_t Cigar::UnitLength(Cigar::Unit const& unit, bool use_top_coord)
{
    return use_top_coord 
        ? unit.op.ref * unit.length
        : unit.op.temp * unit.length;

}


int64_t Cigar::UnitOffset(Cigar::Unit const& unit, bool top_to_bottom)
{
    int multiplier = top_to_bottom ? 1 : -1;
    return (unit.op.ref - unit.op.temp) * unit.length * multiplier;
}



//Take a substring of the CIGAR virtually
Cigar::CIGAR_VEC Cigar::Substring(Cigar::CIGAR_VEC const& cigar, size_t start_position,
                                  size_t end_position, bool use_top_coord)
{
    Cigar::CIGAR_VEC subcigar;
    size_t cigar_pos = 0, prev_cigar_pos = 0;

    if (start_position == end_position)
    {
        return subcigar;
    }

    if (start_position > end_position)
    {
        fprintf(stderr, "Substring: start_position > end_position\n");
        exit(1);
    }

    //count along the top sequence length
    for (size_t c = 0; c != cigar.size(); ++c)
    {
        Cigar::Unit const& unit = cigar[c];
        cigar_pos += Cigar::UnitLength(unit, use_top_coord);
        
        if (cigar_pos > start_position)
        {
            subcigar.push_back(unit);
            if (subcigar.size() == 1)
            {
                subcigar[0].length -= (start_position - prev_cigar_pos);
            }
        }
        
        if (cigar_pos >= end_position)
        {
            subcigar[subcigar.size() - 1].length -= (cigar_pos - end_position);
            break;
        }
        prev_cigar_pos = cigar_pos;
    }
    return subcigar;
}


size_t Cigar::Length(Cigar::CIGAR_VEC const& cigar, bool use_top_coord)
{
    size_t length = 0;
    for (size_t c = 0; c != cigar.size(); ++c)
    {
        length += UnitLength(cigar[c], use_top_coord);
    }
    return length;

}


size_t Cigar::Length(Cigar::CIGAR_VEC::const_iterator cigar_start,
                     Cigar::CIGAR_VEC::const_iterator cigar_end,
                     bool use_top_coord)
{
    size_t length = 0;
    Cigar::CIGAR_VEC::const_iterator unit_iter;
    
    for (unit_iter = cigar_start; unit_iter != cigar_end; ++unit_iter)
    {
        length += UnitLength(*unit_iter, use_top_coord);
    }
    return length;

}


//return the 
size_t Cigar::Overlap(Cigar::CIGAR_VEC::const_iterator cigar_start,
                      Cigar::CIGAR_VEC::const_iterator cigar_end)
{
    Cigar::CIGAR_VEC::const_iterator unit_iter;
    size_t length = 0;
    
    for (unit_iter = cigar_start; unit_iter != cigar_end; ++unit_iter)
    {
        if ((*unit_iter).op.code == Cigar::M)
        {
            length += (*unit_iter).length;
        }
    }
    return length;
}

                      
//returns the displacement (offset) from the starts of one sequence to the other.
//if 'top_start_to_bottom_start' is true, offset will be positive iff top < bottom,
//for example.
int64_t Cigar::LeftOffset(CIGAR_VEC const& cigar, bool top_start_to_bottom_start)
{
    Cigar::CIGAR_VEC::const_iterator unit_iter;

    int64_t sense = top_start_to_bottom_start ? 1 : -1;
    int64_t offset = 0;
    for (unit_iter = cigar.begin(); unit_iter != cigar.end(); ++unit_iter)
    {
        Unit const& u = (*unit_iter);
        if (u.op.code == Cigar::M)
        {
            break;
        }
        else
        {
            offset += (u.op.ref - u.op.temp) * u.length;
        }
    }
    return offset * sense;
}


int64_t Cigar::RightOffset(CIGAR_VEC const& cigar, bool top_end_to_bottom_end)
{
    Cigar::CIGAR_VEC::const_reverse_iterator unit_iter;

    int64_t sense = top_end_to_bottom_end ? -1 : 1;
    int64_t offset = 0;
    for (unit_iter = cigar.rbegin(); unit_iter != cigar.rend(); ++unit_iter)
    {
        Unit const& u = (*unit_iter);
        if (u.op.code == Cigar::M)
        {
            break;
        }
        else
        {
            offset += (u.op.ref - u.op.temp) * u.length;
        }
    }
    return offset * sense;
}


/*
We traverse both CIGARs block by block, producing a third, transitive merge CIGAR based on the
lengths and CIGAR operations of the comparison.

The matrix cells are:

merged CIGAR op, merged CIGAR length
amount-to-advance first, amount-to-advance second

MM  M
MI  I
MD  D
IM  D
II  error
ID  D
DM  I
DI  I
DD  <disappears>

NM  
NI
ND
         2nd CIGAR
   M         I        D

M  M,min     I,2nd    D,min
   min,min   0,2nd    min,min

I  D,1st     error    D,1st
   1st,0     ----     1st,0

D  I,min     I,2nd    *
   min,min   0,2nd    min,min

 */

//sorted by 
std::multiset<std::pair<size_t, size_t> > Cigar::ComputeOffsets(Cigar::CIGAR_VEC const& cigar)
{
    std::multiset<std::pair<size_t, size_t> > cumul_offsets;
    size_t cumul_top_offset = 0;
    size_t cumul_bot_offset = 0;

    for (Cigar::CIGAR_ITER iter = cigar.begin(); iter != cigar.end(); ++iter)
    {
        cumul_top_offset += Cigar::UnitLength(*iter, true);
        cumul_bot_offset += Cigar::UnitLength(*iter, false);
        cumul_offsets.insert(std::pair<size_t, size_t>(cumul_top_offset, cumul_bot_offset));
    }
    return cumul_offsets;
}


Cigar::CIGAR_VEC
Cigar::Expand(std::vector<block_offsets> const& blocks,
              Cigar::CIGAR_VEC const& cigar,
              bool use_N_as_insert)
{
    Cigar::CIGAR_VEC result;

    if (cigar.empty())
    {
        return result;
    }
    
    assert(! blocks.empty());

    Cigar::Op ins = use_N_as_insert ? Ops[Cigar::N] : Ops[Cigar::D];

    //initialize
    std::vector<block_offsets>::const_iterator cur_block = blocks.begin();
    CIGAR_VEC::const_iterator cur_op = cigar.begin();

    size_t op_remain = (*cur_op).length;
    size_t block_remain = (*cur_block).block_length;
    
    while (1)
    {
        if (cur_block == blocks.end())
        {
            //while we are in the loop, there must be at least one
            //block left to define the cigar projection
            fprintf(stderr, "Error: Cigar::Expand: Ran out of blocks to expand this CIGAR\n");
            exit(1);
        }
        
        if (op_remain <= block_remain)
        {
            //Remaining piece of op
            result.push_back(*cur_op);
            block_remain -= (*cur_op).length;
            if (++cur_op == cigar.end())
            {
                break;
            }
            op_remain = (*cur_op).length;
        }
        else
        {
            //CIGAR op is broken up into two pieces
            result.push_back(Cigar::Unit((*cur_op).op, block_remain));
            result.push_back(Cigar::Unit(ins, (*cur_block).jump_length));
            ++cur_block;
            block_remain = (*cur_block).block_length;
        }
    }
    return result;
}


//condenses cigar according to the expansion described in 'blocks'
//unlike with 'Cigar::Expand', here the jumps in 'blocks' are applied
//to the cigar to remove the jumped-over part.
Cigar::CIGAR_VEC
Cigar::Condense(std::vector<block_offsets> const& blocks,
                Cigar::CIGAR_VEC const& cigar)
{
    Cigar::CIGAR_VEC result;

    if (cigar.empty())
    {
        return result;
    }
    
    assert(! blocks.empty());

    //initialize
    std::vector<block_offsets>::const_iterator cur_block = blocks.begin();
    CIGAR_VEC::const_iterator cur_op = cigar.begin();

    size_t op_remain = (*cur_op).length;
    size_t block_remain = (*cur_block).block_length;
    size_t jump_remain = (*cur_block).jump_length;
    
    while (1)
    {
        if (jump_remain > 0)
        {
            //condensing logic.
            //M -> I
            //T -> U
            //else -> nothing.
            size_t min_length = std::min(jump_remain, op_remain);
            
            if ((*cur_op).op.code == Cigar::M)
            {
                result.push_back(Cigar::Unit(Cigar::Ops[Cigar::I], min_length));
            }
            else if ((*cur_op).op.code == Cigar::T)
            {
                result.push_back(Cigar::Unit(Cigar::Ops[Cigar::U], min_length));
            }
            else
            {
                //do nothing.
            }
            jump_remain -= min_length;
            op_remain -= min_length;
        }
        else
        {
            //matching logic
            //merely regurgitate the current op with its remaining length
            size_t min_length = std::min(block_remain, op_remain);
            result.push_back(Cigar::Unit((*cur_op).op, min_length));
            block_remain -= min_length;
            op_remain -= min_length;
        }

        if (block_remain == 0)
        {
            if (++cur_block == blocks.end())
            {
                break;
            }
            jump_remain = (*cur_block).jump_length;
            block_remain = (*cur_block).block_length;
        }
        if (op_remain == 0)
        {
            if (++cur_op == cigar.end())
            {
                break;
            }
            op_remain = (*cur_op).length;
        }
    }
    if (cur_op != cigar.end())
    {
        //while we are in the loop, there must be at least one
        //block left to define the cigar projection
        fprintf(stderr, "Error: Cigar::Condense: Ran out of blocks to condense this CIGAR\n");
        exit(1);
    }

    return result;

}


//returns the cigar relationship CIGAR(first, second) from CIGAR(ref,
//first) and CIGAR(ref, second).  The 'cigar1' may in general be long
//representing an entire transcriptome on one contig, for example.
//It must be indexed, allowing efficient merging of first and second.

/*
Cigar::CIGAR_VEC 
Cigar::TransitiveMerge(Cigar::CIGAR_VEC const& cigar1,
                       std::multiset<std::pair<size_t, size_t> > const& cigar_index1,
                       Cigar::CIGAR_VEC const& cigar2,
                       bool add_padding,
                       bool inserts_are_introns)
{
    Cigar::CIGAR_VEC cigar_trans;

    if (cigar2.empty())
    {
        return cigar_trans;
    }
    
    Cigar::CIGAR_ITER iter1 = cigar1.begin();
    Cigar::CIGAR_ITER end1 = cigar1.end();

    Cigar::CIGAR_ITER iter2 = cigar2.begin();
    Cigar::CIGAR_ITER end2 = cigar2.end();

    Cigar::Op op1 = (*iter1).op;
    size_t remain1 = (*iter1).length;

    Cigar::Op op2 = (*iter2).op;
    size_t remain2 = (*iter2).length;

    //preprocessing step to handle a long D state
    if (op2 == Cigar::D && remain1 < remain2)
    {
        //find the top and bottom length on the transformation of
        //the chunks before the containing chunk
        std::multiset<std::pair<size_t, size_t> >::const_iterator index_iter = 
            cigar_index1.upper_bound(std::pair<size_t, size_t>(remain2, UINT64_MAX));

        assert(index_iter != cigar_index1.begin());

        size_t dist = std::distance(cigar_index1.begin(), index_iter);

        //set the next chunk to be *after 
        std::advance(iter1, dist);

        if (iter1 == end1)
        {
            op1 = Cigar::D;
            remain1 = UINT64_MAX;
        }
        else
        {
            op1 = (*iter1).op;
            remain1 = (*iter1).length;
        }

        //we want to have the cumulative distances up to the point before this
        --index_iter;
        
        //   now index_iter points to the index element that tells us
        //   the part of the transformation falling entirely in the first
        //   test 'D' state.  The net deleted in the transform up to this
        //   point, if padding is enabled, should be added to the
        //   padding.

        if (add_padding)
        {
            size_t pad = (*index_iter).first - (*index_iter).second;
            cigar_trans.push_back(Cigar::Unit(Cigar::P, pad));
        }

        cigar_trans.push_back(Cigar::Unit(Cigar::D, (*index_iter).second));
        remain2 -= (*index_iter).first;
    }
    
    Cigar::Op delete_op = inserts_are_introns ? Cigar::N : Cigar::D;

    size_t minstep;

    //bool cigar1_finished = (iter1 == end1);
    bool cigar2_finished = (iter2 == end2);

    while (! cigar2_finished)
    {

        //if op2 == D, gobble up op1 units until the unit that straddles the
        //end coordinate of op2
        minstep = std::min(remain1, remain2);
        
        switch (op1)
        {
        case Cigar::M:
            switch (op2)
            {
            case Cigar::M:
            case Cigar::D:
            case Cigar::N:
                {
                    cigar_trans.push_back(Cigar::Unit(op2, minstep));
                    remain1 -= minstep;
                    remain2 -= minstep;
                }
                break;
            case Cigar::I:
            case Cigar::S:
            case Cigar::H:
                {
                    cigar_trans.push_back(Cigar::Unit(op2, remain2));
                    remain2 = 0;
                }
                break;
            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;
        case Cigar::I:
            switch (op2)
            {
            case Cigar::M:
            case Cigar::D:
            case Cigar::N:
            case Cigar::S:
            case Cigar::H:
                {
                    cigar_trans.push_back(Cigar::Unit(delete_op, remain1));
                    remain1 = 0;
                }
                break;
            case Cigar::I:
                {
                    fprintf(stderr, "CigarTransitiveMerge: Cannot combine two I states."
                            "Undefined alignment\n");
                    exit(1);
                }
                break;
            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;
        case Cigar::D:
            switch (op2)
            {
            case Cigar::M:
                cigar_trans.push_back(Cigar::Unit(Cigar::I, minstep));
                remain1 -= minstep;
                remain2 -= minstep;
                break;

            case Cigar::I:
                cigar_trans.push_back(Cigar::Unit(Cigar::I, remain2));
                remain2 = 0;
                break;

            case Cigar::D:
            case Cigar::N:
                if (add_padding)
                {
                    cigar_trans.push_back(Cigar::Unit(Cigar::P, minstep));
                }
                remain1 -= minstep;
                remain2 -= minstep;
                break;

            case Cigar::S:
            case Cigar::H:
                {
                    //ignore this operation in the second block
                    remain2 = 0;
                }
                break;

            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;

        default:
            fprintf(stderr, "Error: transitive merge only valid for M, I, and D CIGAR ops.");
            exit(2);
        }
            
        if (remain1 == 0)
        {
            ++iter1;
            if (iter1 != end1)
            {
                remain1 = (*iter1).length;
                op1 = (*iter1).op;
            }
            else
            {
                op1 = Cigar::D;
                remain1 = UINT64_MAX;
            }
        }
        if (remain2 == 0)
        {
            ++iter2;
            cigar2_finished = (iter2 == end2);
            if (! cigar2_finished)
            {
                remain2 = (*iter2).length;
                op2 = (*iter2).op;
            }
        }
        if (cigar2_finished)
        {
            break;
        }
    }

    //at least one of the CIGARs is used up
    // if (cigar1_finished != cigar2_finished)
    // {
    //     //add in the extra.  When a CIGAR is used up, we pretend that it has 
    //     //a long D state at the end.  We know that a D state combines as
    //     //

    //     size_t extra;
    //     Cigar::Op op;
    //     if (cigar1_finished)
    //     {
    //         //this should not happen
    //         assert(false);
    //         // extra = Cigar::UnitLength(Cigar::Unit(op2, remain2), false)
    //         //     + Cigar::Length(iter2 + 1, end2, false);
    //         // op = Cigar::S;
    //         // cigar_trans.push_back(Cigar::Unit(op, extra));
    //     }
    //     else
    //     {
    //         //cigar2 is finished first.  test is shorter.  do nothing.
    //         // extra = Cigar::UnitLength(Cigar::Unit(op1, remain1), false)
    //         //     + Cigar::Length(iter1 + 1, end1, false);
    //         // op = Cigar::D;
    //         // cigar_trans.push_back(Cigar::Unit(op, extra));
    //     }
    // }
    
    //now consolidate consecutive like-states in the CIGAR vector
    return Cigar::Consolidate(cigar_trans);

}
*/


//Used for discovering the agreeing subset of two very similar alignments.
//inputs: cigar1, (CIGAR(ref, first))
//        cigar2, (CIGAR(ref,second))
//output: tcigar: (CIGAR(ref, trimmed_first))
//'M' states that are the intersection of cigar1 and cigar2
//are retained in tcigar, the remaining states are converted to 'S'

 /*
Cigar::CIGAR_VEC 
Cigar::TransitiveTrim(Cigar::CIGAR_VEC const& cigar1,
                      Cigar::CIGAR_VEC const& cigar2,
                      bool add_padding,
                      bool inserts_are_introns)
{
    Cigar::CIGAR_VEC cigar_trans;

    if (cigar2.empty())
    {
        return cigar_trans;
    }
    
    Cigar::CIGAR_ITER iter1 = cigar1.begin();
    Cigar::CIGAR_ITER end1 = cigar1.end();

    Cigar::CIGAR_ITER iter2 = cigar2.begin();
    Cigar::CIGAR_ITER end2 = cigar2.end();

    Cigar::Op op1 = (*iter1).op;
    size_t remain1 = (*iter1).length;

    Cigar::Op op2 = (*iter2).op;
    size_t remain2 = (*iter2).length;

    Cigar::Op delete_op = inserts_are_introns ? Cigar::N : Cigar::D;

    size_t minstep;

    //bool cigar1_finished = (iter1 == end1);
    bool cigar2_finished = (iter2 == end2);

    while (! cigar2_finished)
    {

        //if op2 == D, gobble up op1 units until the unit that straddles the
        //end coordinate of op2
        minstep = std::min(remain1, remain2);
        
        switch (op1)
        {
        case Cigar::M:
            switch (op2)
            {
            case Cigar::M:
            case Cigar::D:
            case Cigar::N:
                {
                    cigar_trans.push_back(Cigar::Unit(op2, minstep));
                    remain1 -= minstep;
                    remain2 -= minstep;
                }
                break;
            case Cigar::I:
            case Cigar::S:
            case Cigar::H:
                {
                    cigar_trans.push_back(Cigar::Unit(op2, remain2));
                    remain2 = 0;
                }
                break;
            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;
        case Cigar::I:
            switch (op2)
            {
            case Cigar::M:
            case Cigar::D:
            case Cigar::N:
            case Cigar::S:
            case Cigar::H:
                {
                    cigar_trans.push_back(Cigar::Unit(delete_op, remain1));
                    remain1 = 0;
                }
                break;
            case Cigar::I:
                {
                    fprintf(stderr, "CigarTransitiveMerge: Cannot combine two I states."
                            "Undefined alignment\n");
                    exit(1);
                }
                break;
            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;
        case Cigar::D:
            switch (op2)
            {
            case Cigar::M:
                cigar_trans.push_back(Cigar::Unit(Cigar::I, minstep));
                remain1 -= minstep;
                remain2 -= minstep;
                break;

            case Cigar::I:
                cigar_trans.push_back(Cigar::Unit(Cigar::I, remain2));
                remain2 = 0;
                break;

            case Cigar::D:
            case Cigar::N:
                if (add_padding)
                {
                    cigar_trans.push_back(Cigar::Unit(Cigar::P, minstep));
                }
                remain1 -= minstep;
                remain2 -= minstep;
                break;

            case Cigar::S:
            case Cigar::H:
                {
                    //ignore this operation in the second block
                    remain2 = 0;
                }
                break;

            case Cigar::P:
            case Cigar::None:
                assert(false);
                break;
            }
            break;

        default:
            fprintf(stderr, "Error: transitive merge only valid for M, I, and D CIGAR ops.");
            exit(2);
        }
            
        if (remain1 == 0)
        {
            ++iter1;
            if (iter1 != end1)
            {
                remain1 = (*iter1).length;
                op1 = (*iter1).op;
            }
            else
            {
                op1 = Cigar::D;
                remain1 = UINT64_MAX;
            }
        }
        if (remain2 == 0)
        {
            ++iter2;
            cigar2_finished = (iter2 == end2);
            if (! cigar2_finished)
            {
                remain2 = (*iter2).length;
                op2 = (*iter2).op;
            }
        }
        if (cigar2_finished)
        {
            break;
        }
    }
//now consolidate consecutive like-states in the CIGAR vector
return Cigar::Consolidate(cigar_trans);
}
*/


  /*
  //how to tell if this is a valid projection?
  size_t Cigar::ProjectCoord(Cigar::CIGAR_VEC const& transformation,
  std::multiset<std::pair<size_t, size_t> > const& trans_index,
  size_t source_coord,
                           bool * projection_applied)
{

    Cigar::CIGAR_VEC source_cigar(1, Cigar::Unit(Cigar::Ops[Cigar::D], source_coord));
    source_cigar.push_back(Cigar::Unit(Cigar::Ops[Cigar::M], 1));

    Cigar::CIGAR_VEC source_proj = 
        Cigar::TransitiveMerge(transformation, trans_index, source_cigar, false, false);

    Cigar::CIGAR_ITER trimmed_left, trimmed_right;
    Cigar::Trim(source_proj, false, &trimmed_left, &trimmed_right);

    size_t overlap = Cigar::Overlap(trimmed_left, trimmed_right);
    
    *projection_applied = (overlap > 0);

    return Cigar::LeftOffset(source_proj, true);
}
*/


Cigar::CIGAR_VEC Cigar::Consolidate(Cigar::CIGAR_VEC const& source)
{
    Cigar::CIGAR_VEC condensed;
    if (source.empty())
    {
        return condensed;
    }
    //condensed.reserve(source.size());
    condensed.push_back(source[0]);

    Cigar::Op prev_op = source[0].op;

    for (size_t s = 1; s != source.size(); ++s)
    {
        if (source[s].op.code == prev_op.code)
        {
            (*condensed.rbegin()).length += source[s].length;
        }
        else
        {
            condensed.push_back(source[s]);
        }
        prev_op = source[s].op;
    }
    return condensed;
}



void Cigar::Trim(Cigar::CIGAR_VEC const& source, 
            bool retain_top_extent,
            CIGAR_ITER * trimmed_start,
            CIGAR_ITER * trimmed_end)
{

    Cigar::CIGAR_VEC::const_iterator left = source.begin();
    Cigar::CIGAR_VEC::const_reverse_iterator right = source.rbegin();
    

    while (left != source.end()
           && Cigar::UnitLength(*left, retain_top_extent) == 0)
    {
        ++left;
    }
    Cigar::CIGAR_VEC::const_reverse_iterator right_end(left);
    
    while (right != right_end
           && Cigar::UnitLength(*right, retain_top_extent) == 0)
    {
        ++right;
    }
    *trimmed_start = left;
    *trimmed_end = right.base();
}


Cigar::CIGAR_VEC Cigar::Trim(Cigar::CIGAR_VEC const& source, 
                             bool retain_top_extent)
{
    Cigar::CIGAR_VEC::const_iterator left, right;
    Trim(source, retain_top_extent, &left, &right);
    return Cigar::CIGAR_VEC(left, right);
}


/*
Cigar::CIGAR_VEC Cigar::Invert(Cigar::CIGAR_VEC const& source, size_t position)
{
    Cigar::CIGAR_VEC invert;
    Cigar::Op inverse_ops[] = { Cigar::M, Cigar::D, Cigar::I, 
                                Cigar::I, Cigar::None, Cigar::None, 
                                Cigar::P };

    for (Cigar::CIGAR_ITER sit = source.begin();
         sit != source.end(); ++sit)
    {
        Cigar::Op inv = inverse_ops[(*sit).op];
        if (inv == Cigar::None)
        {
            fprintf(stderr, "Couldn't invert CIGAR string containing 'S' or 'H'\n");
            exit(1);
        }
        invert.push_back(Cigar::Unit(inv, (*sit).length));
    }
    return invert;
}
*/


//count the number of correctly aligned positions, represented by a
//merged cigar.
size_t Cigar::CountAlignedPositions(Cigar::CIGAR_VEC const& merged, 
                                    size_t * total_bases_top, 
                                    size_t * total_bases_bottom)
{
    size_t guide_read_pos = 0;
    size_t test_read_pos = 0;

    size_t num_correct_bases = 0;
    //once the merge is complete, there will (in general) be an offset between the guide
    //and test start (represented by an 'I' or 'D' state).  
    for (Cigar::CIGAR_ITER mi = merged.begin(); mi != merged.end(); ++mi)
    {
                
        Cigar::Unit const& unit = *mi;

        switch(unit.op.code)
        {
        case Cigar::M:
            if (test_read_pos == guide_read_pos)
            {
                num_correct_bases += unit.length;
            }                    
            else
            {
            }
            guide_read_pos += unit.length;
            test_read_pos += unit.length;

            break;

        case Cigar::I:
        case Cigar::U:
            //present in test but not in guide.  test is erroneously placed
            test_read_pos += unit.length;
            break;
            
        case Cigar::D:
        case Cigar::N:
            guide_read_pos += unit.length;
            break;
                    
        case Cigar::H:
            assert(false);
            break;

        case Cigar::P:
            break;

        case Cigar::S:
            //this shouldn't happen since we're dealing with a merged CIGAR
            assert(false);
            break;

        case Cigar::T:
            guide_read_pos += unit.length;
            test_read_pos += unit.length;
            break;

        }
    }
    *total_bases_top = guide_read_pos;
    *total_bases_bottom = test_read_pos;

    return num_correct_bases;
}

//for the projection to be well-defined, the coordinate should
//fall in an M block.

//basic idea: traverse the CIGAR, updating the source and target positions.
//as soon as you reach the coordinate, either:
//1. you exceed the coordinate, and thus, the current block must be M
//2. you hit the coordinate exactly.  either the previous or current block must be M

/*
int64_t Cigar::ProjectCoord(Cigar::CIGAR_ITER projection_start,
                            Cigar::CIGAR_ITER projection_end,
                            int64_t coordinate,
                            bool top_to_bottom,
                            bool * is_missing_projection)
{
    
    int64_t source_position = 0;
    int64_t target_position = 0;
    Cigar::CIGAR_ITER op;
    for (op = projection_start; op != projection_end; ++op)
    {
        source_position += Cigar::UnitLength(*op, top_to_bottom);
        target_position += Cigar::UnitLength(*op, ! top_to_bottom);

        if (coordinate <= source_position)
        {
            //we've covered enough of the CIGAR to encounter the coordinate
            if (coordinate < source_position)
            {
                if ((*op).op != Cigar::M)
                {
                    * is_missing_projection = true;
                    //this projection is in a 'meta-intron'.
                    //report the projected position as if it were at the beginning of the
                    //last valid target position
                    return target_position;
                }
            }
            else
            {
                //coordinate == source_position, we should be on a boundary
            }
            break;
        }
    }
    //we're in a block such that start < coordinate <= end
    //this must be a M block, otherwise this is undefined behavior
    *is_missing_projection = false;
    return coordinate + target_position - source_position;
}
*/


