#ifndef _SEQ_PROJECTION_H
#define _SEQ_PROJECTION_H

#include "sam_helper.h"

#include <vector>

struct block_offsets
{
    size_t jump_length;
    size_t block_length;
    block_offsets(size_t _j, size_t _b);
};


std::vector<block_offsets> InitFromMDString(char const* md_string);

class SequenceProjection
{
 public:
    std::string species;
    std::string source_dna;
    std::string target_dna;
    bool same_strand;
    std::vector<block_offsets> transformation;

    size_t total_block_length;

    /* Cigar::CIGAR_VEC cigar; */
    /* Cigar::CIGAR_INDEX cigar_index; */

    SequenceProjection(char const* _species,
                       char const* _source_dna,
                       char const* _target_dna,
                       char _strand_char,
                       std::vector<block_offsets> const& blocks);

    bool operator<(SequenceProjection const& b) const;
};


bool ApplyProjectionToSAM(SequenceProjection const& projection,
                          char const* alignment_space,
                          SamLine * samline,
                          bool inserts_are_introns,
                          bool add_cufflinks_xs_tag);

#endif // _SEQ_PROJECTION_H
