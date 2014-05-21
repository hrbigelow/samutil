#ifndef _SAM_AUX_H
#define _SAM_AUX_H

#include <unordered_map>
#include "sam_score_aux.h"
#include "align_eval_aux.h"


// typedef std::set<SequenceProjection, less_seq_projection>::const_iterator SEQ_PROJ_ITER;

typedef std::unordered_map<char const*, SEQ_PROJ_ITER, to_integer, eqstr> PROJ_MAP;

typedef PROJ_MAP::const_iterator NP_ITER;


/* SEQ_PROJ_ITER  */
/* non_overlapping_range(SEQ_PROJ_ITER start, SEQ_PROJ_ITER set_end); */


struct pdp_result
{
    pdp_result();
    std::vector<char> * lines;
    size_t num_records_retained;
    size_t num_records_discarded;
};

// inputs: SamLine * [beg, end) range, SamBuffer *
// outputs: char * projected_deduped_samlines
// allocates output.
struct project_dedup_print
{
    SamOrder const* genome_sam_order;
    bool inserts_are_introns;
    bool retain_unsequenced_projection;
    size_t average_rsam_line_length;

    project_dedup_print(SamOrder const* _genome_order,
                        bool _inserts_are_introns,
                        bool _retain_unseq,
                        size_t _ave_len);

    pdp_result operator()(std::pair<SAMIT, SAMIT> const& range);

};


// inputs: rSAM-formatted SamLine *, INDEX_ITER
// outputs: SAM-formatted allocated SAM record string
struct rsam_to_sam_binary
{
    char const* seq_buffer; // generated from 'samutil seqindex'
    size_t data_buffer_offset; // the start_offset corresponding to the beginning of this seq_buffer
    SamFilter const* sam_filter;

    rsam_to_sam_binary(char const* _seq_buffer, size_t _dbo, SamFilter const* _sf);
    char * operator()(SamLine * samline, INDEX_ITER li_iter);

};


// inputs: rSAM-formatted SamLine *, INDEX_ITER
// outputs: SAM-formatted allocated SAM record string
struct unmapped_rsam_to_fastq
{
    char const* seq_buffer; // generated from 'samutil seqindex'
    size_t data_buffer_offset; // the start_offset corresponding to the beginning of this seq_buffer

    unmapped_rsam_to_fastq(char const* _seq_buffer, size_t _dbo);
    char const* operator()(SamLine * samline, INDEX_ITER li_iter);

};


struct delete_samline
{
    delete_samline();
    void operator()(SamLine * samline);
};


// assume a tab-separated line of [id seq1 qual1 [seq2 qual2 [...]]]
void print_fqd_as_fastq(char const* fqd_line, FILE * fq_fh);


size_t average_line_length(std::vector<SamLine *>::const_iterator beg,
                           std::vector<SamLine *>::const_iterator end,
                           size_t nsample);


#endif // _SAM_AUX_H
