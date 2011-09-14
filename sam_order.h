#ifndef _SAM_ORDER_H
#define _SAM_ORDER_H

#include "sam_helper.h"
#include "gtf.h"

#include <map>
#include <set>
#include <string>

enum SAM_ORDER
    {
        SAM_RID,
        SAM_POSITION_RID,
        SAM_RID_POSITION,
        SAM_POSITION,
        SAM_FID_POSITION
    };

/* enum SAM_INDEX_TYPE */
/*     { */
/*         SAM_POS_MIN_ALIGN_GUIDE, */
/*         SAM_POS_ALIGN, */
/*         SAM_READ_ID, */
/*         SAM_READ_ID_ILLUMINA */
/*     }; */


class SamOrder
{
 public:
    SamOrder(SAM_ORDER _so, char const* _sort_type);
    SamOrder();
    ~SamOrder();
    SAM_QNAME_FORMAT InitFromFile(FILE * sam_fh);
    SAM_QNAME_FORMAT InitFromFastqFile(char const* fq_file);
    SAM_QNAME_FORMAT InitFromID(char const* id);

    void InitFromChoice(SAM_QNAME_FORMAT qname_format);

    SAM_ORDER order;
    /* SAM_INDEX_TYPE index_type; */

    char sort_type[50];

    std::map<std::string, size_t> contig_lengths;
    CONTIG_OFFSETS contig_offsets;

    bool (SamOrder::* less)(SamLine const& a, SamLine const& b) const;
    bool (SamOrder::* equal)(SamLine const& a, SamLine const& b) const;
    size_t (* parse_fragment_id)(char const* qname);

    size_t (SamOrder::* sam_index)(char const* samline) const;
    
    
    void AddHeaderContigStats(FILE * sam_fh);

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

    size_t samline_position_min_align_guide(char const* samline) const;
    size_t samline_position_align(char const* samline) const;
    size_t samline_read_id(char const* samline) const;
    size_t samline_fragment_id(char const* samline) const;

};


void GenerateProjectionHeader(CONTIG_OFFSETS const& genome_contig_order,
                              std::set<SequenceProjection> const& tx_to_genome,
                              FILE * out_sam_header_fh);

size_t parse_fragment_id_numeric(char const* qname);

size_t parse_fragment_id_illumina(char const* qname);

size_t parse_fragment_id_casava_1_8(char const* qname);

size_t parse_fragment_id_zero(char const* /* qname */);

#endif // _SAM_ORDER_H
