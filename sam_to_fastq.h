/*
  Utilities for parsing a SAM file, ordered by FRAGMENT, into a set of
  fastq files.

  Assumes input is sorted by FRAGMENT (which is unspecified regarding read pair number)
  
  User provides number of reads in fragment. (According to SAM spec,
  the number of 'segments' in the 'template'). Probably this will only support '1' or '2'.

  Outputs R + 1 fastq(.gz) files, one for each read number, plus one for 'orphan' reads.

  Reads / Quals are reverse-complemented relative to the SAM output as necessary to restore the
  original fastq data.

  Expects a std::vector<char *> of zero-terminated sam_lines, plus buffers for allocated memory

 */

#include "sam_order.h"

// a view of a SAM line entry as Fastq
struct SAMFastqView
{
    char const* qname;
    size_t qname_len;
    size_t fragment_id;
    SamFlag flag;
    char const* seq;
    char const* qual;
    size_t seqlen;
    bool do_print;
    bool is_orphan;

    SAMFastqView();
    SAMFastqView(char const* _qname, size_t _qname_len,
                 size_t _fragment_id, size_t _flag_raw, 
                 char const* _seq, char const* _qual, 
                 size_t _seqlen, bool _do_print, bool _is_orphan);
    char * print(char *outbuf);
};


// very similar to 'partial_index_aux' from align_eval_aux.h
struct SAMFastqView_aux
{
    SamOrder * sam_order;
    SAMFastqView_aux(SamOrder * _sam_order);
    SAMFastqView operator()(char * samline_string);
};


// convert samlines to fastq, filtering out duplicates
// and properly classifying as fastq1, fastq2, or orphan
// if fastq2_buffer_out is NULL, assume we have single-end reads
// otherwise assume paired-end reads
// return the number of samlines NOT processed.
size_t convert_chunk_sam_to_fastq(std::vector<char *> & samlines,
                                  SamOrder const* sam_order,
                                  char * chunk_buffer_in,
                                  char * fastq1_buffer_out,
                                  char * fastq2_buffer_out,
                                  char * orphan_buffer_out,
                                  bool is_last_chunk);

void init_fastq_view(std::vector<char *> & sam_lines,
                     SAMFastqView * fastq_view);


size_t find_last_fragment_bound(SAMFastqView * fastq_view, size_t S, bool is_last_chunk);

void set_flags_fastq_view_paired(std::vector<char *> & sam_lines,
                                 SAMFastqView * fastq_view,
                                 bool is_last_chunk);


void set_flags_fastq_view_single(std::vector<char *> & sam_lines,
                                 SAMFastqView * fastq_view,
                                 bool is_last_chunk);


void write_fastq_view_paired(SAMFastqView * fastq_view,
                             size_t vS,
                             char * fastq1_buffer_out,
                             char * fastq2_buffer_out,
                             char * orphan_buffer_out,
                             size_t * nbytes_written);


void write_fastq_view_single(SAMFastqView * fastq_view,
                             size_t vS,
                             char * fastq_buffer_out,
                             size_t * nbytes_written);
