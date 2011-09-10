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



//apply a partial projection to rSAM line
//this is a helper function for ApplyProjectionToSAM
//returns true if projection successfully applied
bool apply_projection_aux(SequenceProjection const& projection,
                          SamLine * samline,
                          bool inserts_are_introns)
{
    if (samline->flag.this_fragment_unmapped)
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
        char bufstring[1024];
        Cigar::ToString(trimmed_left, trimmed_right, bufstring);
        assert(samline->extra == NULL);
        assert(strlen(bufstring) < 1024);

        size_t new_field_size = 
            strlen(bufstring) 
            + projection.target_dna.size()
            + 3;

        samline->extra = new char[new_field_size];
        samline->cigar = samline->extra;
        samline->rname = samline->extra + strlen(bufstring) + 1;

        strcpy(samline->cigar, bufstring);
        samline->pos = expanded_start_offset;
        strcpy(samline->rname, projection.target_dna.c_str());
    }
    return overlap > 0;

}



bool ApplyProjectionToSAM(SequenceProjection const& projection,
                          char const* alignment_space,
                          SamLine * samline,
                          bool inserts_are_introns,
                          bool add_cufflinks_xs_tag)
{
    bool projected = apply_projection_aux(projection, samline, inserts_are_introns);

    samline->add_tag(AlignSpaceTag, AlignSpaceType, alignment_space);
    
    if (add_cufflinks_xs_tag)
    {
        char const* sense = projection.same_strand ? "+" : "-";
        
        samline->add_tag("XS", 'A', sense);
    }
    return projected;
}
