#ifndef _SAM_ORDER_H
#define _SAM_ORDER_H

#include "sam_helper.h"
#include "seq_projection.h"
#include "gtf.h"

#include <map>
#include <set>
#include <string>
#include <omp.h>

enum SAM_ORDER
    {
        SAM_RID,
        SAM_POSITION_RID,
        SAM_RID_POSITION,
        SAM_POSITION,
        SAM_FID_POSITION
    };

#define QNAME_FORMAT_ERROR SIZE_MAX


class SamOrder
{
 public:
    SamOrder(SAM_ORDER _so, char const* _sort_type);
    SamOrder();
    ~SamOrder();
    SAM_QNAME_FORMAT InitFromFile(FILE * sam_fh);
    SAM_QNAME_FORMAT InitFromID(char const* id);

    void InitFromChoice(SAM_QNAME_FORMAT qname_format);
    void InitProjection(char const* gtf_file);

    bool Initialized() const;

    SAM_ORDER order;
    /* SAM_INDEX_TYPE index_type; */

    char sort_type[50];

    std::map<std::string, size_t> contig_lengths;
    CONTIG_OFFSETS contig_offsets;

    PROJECTIONS projections;

    // initialized in constructor
    bool (SamOrder::* less)(SamLine const& a, SamLine const& b) const;
    bool (SamOrder::* equal)(SamLine const& a, SamLine const& b) const;
    size_t (SamOrder::* sam_index)(char const* samline);

    // initialized in SamOrder::InitFromChoice
    size_t (SamOrder::* parse_fragment_id)(char const* qname);
    
    void AddHeaderContigStats(char const* header);

    //gives the flattened position of this alignment according to the
    //concatenated set of contigs in the order specified by contig_offsets.
    //sets contig_iter to point to the found contig.
    size_t flattened_position(SamLine const* a, CONTIG_OFFSETS::const_iterator * contig_iter) const;
    size_t flattened_position_mate(SamLine const* b, CONTIG_OFFSETS::const_iterator * contig_iter) const;

    //uses fields [rname, pos, strand, cigar]
    bool less_position(SamLine const& a, SamLine const& b) const;
    bool equal_position(SamLine const& a, SamLine const& b) const;

    //[minimum of flattened position, query_strand, cigar, mrnm, mpos]
    bool less_fposition(SamLine const& a, SamLine const& b) const;

    //uses fields [rname, first_read_in_pair()]
    bool less_rid(SamLine const& a, SamLine const& b) const;
    bool equal_rid(SamLine const& a, SamLine const& b) const;

    //sort order [qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
    bool less_rid_position(SamLine const& a, SamLine const& b) const;
    bool equal_rid_position(SamLine const& a, SamLine const& b) const;

    // sort order [rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
    bool less_position_rid(SamLine const& a, SamLine const& b) const;
    bool equal_position_rid(SamLine const& a, SamLine const& b) const;

    // sort order [rname]
    bool less_fragment_id(SamLine const& a, SamLine const& b) const;
    bool equal_fragment_id(SamLine const& a, SamLine const& b) const;

    // sort order [rname, pos, strand, cigar]
    bool less_fragment_id_position(SamLine const& a, SamLine const& b) const;
    bool equal_fragment_id_position(SamLine const& a, SamLine const& b) const;

    size_t samline_position_min_align_guide(char const* samline);
    size_t samline_position_align(char const* samline);
    size_t samline_projected_position_align(char const* samline);
    size_t samline_fragment_id(char const* samline);

    /* size_t samline_read_id(char const* samline) const; */

    uint flowcell_hash_value(char const* flowcell);

    // maintains a record of all flowcells seen, and assigns integers to them.
    omp_lock_t flowcell_hash_write_lock;
    std::map<char const*, uint, less_char_ptr> flowcell_hash;

    size_t parse_fragment_id_numeric(char const* qname);
    size_t parse_fragment_id_illumina(char const* qname);
    size_t parse_fragment_id_casava_1_8(char const* qname);
    size_t parse_fragment_id_zero(char const* /* qname */);


};


SAM_QNAME_FORMAT QNAMEFormat(char const* sam_dataline);

struct less_seq_projection
{
    SamOrder const* sam_order;

    less_seq_projection(SamOrder const* _so);

    // sort order by target start position
    bool operator()(SequenceProjection const& a,
                    SequenceProjection const& b) const;
};


size_t flattened_position_aux(char const* contig,
                              size_t position, 
                              CONTIG_OFFSETS const& contig_offsets,
                              CONTIG_OFFSETS::const_iterator * contig_iter);


void GenerateProjectionHeader(CONTIG_OFFSETS const& genome_contig_order,
                              std::set<SequenceProjection> const& tx_to_genome,
                              FILE * out_sam_header_fh);

#endif // _SAM_ORDER_H
