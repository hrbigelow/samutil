#ifndef _SAM_STATS_AUX_H
#define _SAM_STATS_AUX_H

#include "cigar_ops.h"

#include <map>
#include <string>


typedef std::map<std::string, size_t> META_CONTIG;

size_t find_max_consensus_quality(char const* pileup_file);

void samfile_limits(char const* sam_file, 
                    size_t * max_read_length,
                    size_t * max_base_quality,
                    size_t * num_distinct_align_scores);


Cigar::CIGAR_VEC
fill_buffer(FILE * pileup_fh,
            std::map<std::string, size_t> const& contig_offsets,
            Cigar::CIGAR_VEC const& prefix_buf_to_meta_cigar,
            size_t ones_offset,
            size_t fill_size,
            size_t max_consensus_quality,
            int * consensus_score_buf,
            char * consensus_base_buf);


Cigar::CIGAR_VEC 
advance_consensus_buffers(FILE * pileup_fh,
                          std::map<std::string, size_t> const& contig_offsets,
                          Cigar::CIGAR_VEC const& prev_full_cigar,
                          size_t leap_size,
                          size_t ones_offset,
                          size_t max_consensus_quality,
                          int * consensus_score_buf,
                          char * consensus_base_buf);

#endif // _SAM_STATS_AUX_H
