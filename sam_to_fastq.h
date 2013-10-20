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


// a view of a SAM line entry as Fastq
struct SAMFastqView
{
    char const* qname;
    size_t qname_len;
    size_t fragment_id;
    size_t flag_raw;
    char const* seq;
    char const* qual;
    size_t seqlen;
    bool do_print;
    bool is_orphan;

    SAMFastqView();
    SAMFastqView(char const* _qname, size_t _qname_len,
                 size_t _fragment_id, size_t _flag_raw, 
                 char const* _seq, char const* _qual, 
                 size_t _seqlen);
    void print(char *outbuf);
};


// very similar to 'partial_index_aux' from align_eval_aux.h
struct SAMFastqView_aux
{
    SAMOrder const* sam_order;
    SAMFastqView_aux(SamOrder const* _sam_order);
    SAMFastqView operator()(char * samline_string);
};


// convert samlines to fastq, filtering out duplicates
// and properly classifying as fastq1, fastq2, or orphan
// if fastq2_buffer_out is NULL, assume we have single-end reads
// otherwise assume paired-end reads
// return the number of samlines NOT processed.
size_t convert_chunk_sam_to_fastq(std::vector<char *> & samlines,
                                  SAMOrder const* sam_order,
                                  char * chunk_buffer_in,
                                  char * fastq1_buffer_out,
                                  char * fastq2_buffer_out,
                                  char * orphan_buffer_out,
                                  bool is_last_chunk);
