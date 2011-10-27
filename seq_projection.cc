#include "seq_projection.h"

#include <algorithm>

#include "cigar_ops.h"


block_offsets::block_offsets(size_t _j, size_t _b) : jump_length(_j), block_length(_b) { }


std::vector<block_offsets> InitFromMDString(char const* md_string)
{
    std::vector<block_offsets> blocks;
    size_t jump_length = 0;

    char opname;
    int64_t op_length;

    char const* chunk = md_string;
    while (*chunk != '\0')
    {
        size_t num_fields = sscanf(chunk, "%zi%c", &op_length, &opname);
        assert(num_fields == 2);

        chunk = strchr(chunk, opname) + 1;
        assert(chunk != NULL);

        switch(opname)
        {
        case 'M':
            blocks.push_back(block_offsets(jump_length, op_length));
            jump_length = 0;
            break;
        case 'D':
        case 'N':
            //should really be 'jump_length = op_length'.  But, if we
            //encounter consecutive 'D' states, this will handle them
            //correctly.
            jump_length += op_length;
            break;
        default:
            fprintf(stderr, "Error: InitFromMDString: md_string %s "
                    "should be a CIGAR with only M and D operations\n", md_string);
            exit(1);
        }
    }
    return blocks;
}


SequenceProjection::SequenceProjection(char const* _species,
                                       char const* _source_dna,
                                       char const* _target_dna,
                                       char _strand_char,
                                       std::vector<block_offsets> const& _blocks)
{
    species = std::string(_species);
    source_dna = std::string(_source_dna);
    target_dna = std::string(_target_dna);
    same_strand = _strand_char == '+';
    transformation = _blocks;
    total_block_length = 0;
}



SequenceProjection::SequenceProjection(SequenceProjection const& sp) :
    species(sp.species),
    source_dna(sp.source_dna),
    target_dna(sp.target_dna),
    same_strand(sp.same_strand),
    transformation(sp.transformation),
    total_block_length(sp.total_block_length) { }



//order by coordinates of 
bool SequenceProjection::operator<(SequenceProjection const& b) const
{
    
    int64_t a_left_offset = this->transformation.empty() ? 0 : this->transformation[0].jump_length;
    int64_t b_left_offset = b.transformation.empty() ? 0 : b.transformation[0].jump_length;

    // int64_t a_right_neg_offset = - Cigar::RightOffset(this->cigar, true);
    // int64_t b_right_neg_offset = - Cigar::RightOffset(b.cigar, true);

    return this->source_dna < b.source_dna ||
        (this->source_dna == b.source_dna &&
         (a_left_offset < b_left_offset ||
          (a_left_offset == b_left_offset &&
           // (a_right_neg_offset < b_right_neg_offset ||
           //  (a_right_neg_offset == b_right_neg_offset &&
             (this->target_dna < b.target_dna ||
              (this->target_dna == b.target_dna &&
               (this->same_strand < b.same_strand))))));
}


size_t SequenceProjection::target_start_pos() const
{
    assert(! this->transformation.empty());
    return this->transformation[0].jump_length;
}

size_t SequenceProjection::target_end_pos() const
{
    assert(! this->transformation.empty());
    size_t end = 0;
    std::vector<block_offsets>::const_iterator bit;
    for (bit = this->transformation.begin();
         bit != this->transformation.end(); ++bit)
    {
        end += (*bit).jump_length + (*bit).block_length;
    }
    return end;
}


// determine the expanded start position on a chromosome given the
// transcript alignment start position, the CIGAR, and the projection,
// taking strand orientation into account.
size_t ExpandedStartPos(SequenceProjection const& projection,
                        size_t local_start_pos,
                        char const* cigar)
{
    
    size_t start_pos = projection.same_strand 
        ? local_start_pos
        : (projection.total_block_length 
           - Cigar::Length(Cigar::FromString(cigar, 0), true)
           - local_start_pos);

    return ExpandedPos(projection, start_pos);
}


// determine the expanded position relative to a given projection
size_t ExpandedPos(SequenceProjection const& projection,
                   int64_t pos)
{
    if (pos > static_cast<int64_t>(projection.total_block_length))
    {
        fprintf(stderr, "Error: ExpandedPos: position %zi greater than total\n"
                "block length %zu\n", pos, projection.total_block_length);
        exit(1);
    }
    
    size_t ex_pos = pos;
    std::vector<block_offsets>::const_iterator iter = 
        projection.transformation.begin();

    while (iter != projection.transformation.end() && pos >= 0)
    {
        pos -= (*iter).block_length;
        ex_pos += (*iter).jump_length;
        ++iter;
    }
    return ex_pos;
}




//apply a partial projection to rSAM line
//this is a helper function for ApplyProjectionToSAM
//returns true if projection successfully applied
bool ApplySequenceProjection(SequenceProjection const& projection,
                             SamLine * samline,
                             bool inserts_are_introns)
{
    if (! samline->flag.all_fragments_mapped)
    {
        return false;
    }

    Cigar::CIGAR_VEC source_to_read = Cigar::FromString(samline->cigar, 0);

    size_t start_offset = projection.same_strand 
        ? samline->pos
        : (projection.total_block_length 
           - Cigar::Length(source_to_read, true)
           - samline->pos);

    if (! projection.same_strand)
    {
        //reverse alignment.  right_offset is the distance from end of
        //read alignment to end of source.
        std::reverse(source_to_read.begin(), source_to_read.end());

        // flip the template-to-reference, since the reference is
        //changing from transcriptome to genome and they are in
        //opposite directions.
        samline->flag.template_layout ^= 1; 
    }

    size_t expanded_start_offset;

    Cigar::CIGAR_VEC target_to_read =
        Cigar::Expand(projection.transformation, source_to_read, 
                      start_offset, & expanded_start_offset,
                      inserts_are_introns);

    Cigar::CIGAR_ITER trimmed_left, trimmed_right;
    Cigar::Trim(target_to_read, false, &trimmed_left, &trimmed_right);

    size_t overlap = Cigar::Overlap(trimmed_left, trimmed_right);

    if (overlap > 0)
    {
        char cigar_new[1024];
        Cigar::ToString(trimmed_left, trimmed_right, cigar_new);

        char const* rname_new = projection.target_dna.c_str();
        if (strlen(samline->rname) < strlen(rname_new))
        {
            delete samline->rname;
            samline->rname = new char[strlen(rname_new) + 1];
        }
        strcpy(samline->rname, rname_new);

        if (strlen(samline->cigar) < strlen(cigar_new))
        {
            delete samline->cigar;
            samline->cigar = new char[strlen(cigar_new) + 1];
        }
        strcpy(samline->cigar, cigar_new);

        samline->pos = expanded_start_offset;
    }
    return overlap > 0;

}
