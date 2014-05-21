#ifndef _SEQ_PROJECTION_H
#define _SEQ_PROJECTION_H

#include <vector>
#include <unordered_map>
#include <string>
#include <cstddef>

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

    size_t target_start_pos() const;
    size_t target_end_pos() const;

    SequenceProjection(char const* _species,
                       char const* _source_dna,
                       char const* _target_dna,
                       char _strand_char,
                       std::vector<block_offsets> const& blocks);

    SequenceProjection(SequenceProjection const& sp);

    SequenceProjection();

    bool operator<(SequenceProjection const& b) const;
};


size_t ExpandedStartPos(SequenceProjection const& projection,
                        size_t local_start_pos,
                        char const* cigar);


size_t ExpandedPos(SequenceProjection const& projection, int64_t pos);


#endif // _SEQ_PROJECTION_H
