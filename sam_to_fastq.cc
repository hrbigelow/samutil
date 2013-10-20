#include "sam_helper.h"
#include "sam_to_fastq.h"

// do not initialize anything
SAMFastqView::SAMFastqView() { }

// SAMFastqView::SAMFastqView()
//     : qname(NULL), qname_len(0), 
//       fragment_id(0), flag_raw(0), 
//       seq(NULL), qual(NULL), seqlen(0) { }

SAMFastqView::SAMFastqView(char const* _qname, size_t _qname_len, 
                           size_t _fragment_id, size_t _flag_raw, 
                           char const* _seq, char const* _qual, size_t _seqlen) :
    qname(_qname),
    qname_len(_qname_len),
    fragment_id(_fragment_id),
    flag_raw(_flag_raw),
    seq(_seq),
    qual(_qual),
    seqlen(_seqlen)
{ }


// print out a fastq entry.  assume outbuf has enough space
// returns pointer to next insertion point
char * SAMFastqView::print(char * outbuf)
{
    char * w = outbuf;
    *w = '@';
    ++w;
    strncpy(w, this->qname, this->qname_len);
    w += this->qname_len;
    *w = '\n';
    ++w;
    strncpy(w, this->seq, this->seqlen);
    if (this->flag.this_fragment_on_neg_strand)
    {
        reverse_comp(w, w + this->seqlen);
    }
    w += this->seqlen;
    strncpy(w, "\n+\n", 3);
    w += 3;
    strncpy(w, this->qual, this->seqlen);
    if (this->flag.this_fragment_on_neg_strand)
    {
        std::reverse(w, w + this->seqlen);
    }
    w += this->seqlen;
    *w = '\n';
    ++w;
    return w;
}


SAMFastqView_aux::SAMFastqView_aux(SAMOrder const* _sam_order) : sam_order(_sam_order) { }


// Expect a Spec-conformant samline_string, not necessarily zero-terminated
SAMFastqView_aux::operator()(char * samline_string)
{
    size_t fragment_id;
    size_t flag;
    char * qname;
    size_t qname_len;
    char const* seq;
    char const* seq_end;
    char const* qual;
    int seqstart, seqend;
    size_t seqlen;

    // (taken from sam_helper.cc)
    // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
    int num_fields =
        sscanf(samline_string, 
               "%[^\t]%n\t" "%zu\t" 
               "%*[^\t]\t" "%*zu\t" "%*zu\t" "%*s\t" "%*s\t" "%*zu\t" "%*i\t"
               "%n%s%n\t" "%s\t",
               &qname, &qname_len, &flag_raw, &seqstart, &seq, &seqend, &qual);

    if (num_fields != 4)
    {
        fprintf(stderr, "Error: SAMFastqView_aux: couldn't parse samline string into a Fastq view\n");
        exit(35);
    }

    fragment_id = (this->sam_order->*(this->sam_order->sam_index))(qname);

    if (fragment_id == QNAME_FORMAT_ERROR)
    {
        fprintf(stderr, "Error: SAMFastqView_aux: invalid qname (not Illumina format): %s\n", qname)
    }
    seqlen = seqend - seqstart;
    return SAMFastqView(qname, qname_len, fragment_id, flag, seq, qual, seqlen);
}



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
                                  bool is_last_chunk)
{
    bool paired_end_reads = (fastq2_buffer_out != NULL);

    timespec time_begin, time_end;
    
    clock_gettime(CLOCK_REALTIME, &time_begin);
    fprintf(stderr, "Converting to fastq...");
    fflush(stderr);

    // 1. build the index (in parallel)
    size_t S = samlines.size();

    SAMFastqView view = new SAMFastqView[S];
    SAMFastqView * vit;

    SAMFastqView_aux s2v_aux(sam_order);
    __gnu_parallel::transform(samlines.begin(), samlines.end(), view, s2v_aux);

    // 2. traverse the index in a serial loop, dealing and revcomp'ing.
    // find the bound at which we can guarantee complete entries
    size_t last_fragment_id = view[S-1].fragment_id;
    size_t s_last = S;
    if (! is_last_chunk)
    {
        do
        {
            --s_last;
        }
        while (view[s_last].fragment_id == last_fragment_id);
    }
        
    // First pass.  Set flags on 
    if (paired_end_reads)
    {
        size_t current_fragment_id = SIZE_MAX;
        bool read_seen[2] = { false, false };
        size_t marked_read_index[2];
        bool is_orphan;
        size_t read_index; // 0 or 1 depending on whether first in template or not.

        for (vit = view; vit != view + s_last; ++vit)
        {
            read_index = (*vit).flag.first_fragment_in_template ? 0 : 1;
            if (current_fragment_id != (*vit).fragment_id)
            {
                is_orphan = ! (read_seen[0] && read_seen[1]);
                
                // the is_orphan field must wait until seeing all reads of this fragment
                view[marked_read_index[0]].is_orphan = is_orphan;
                view[marked_read_index[1]].is_orphan = is_orphan;

                // reset the state variables for this new fragment
                read_seen[0] = false;
                read_seen[1] = false;
                current_fragment_id = (*vit).fragment_id;
            }
            else
            {
                // the do_print field can be initialized as soon as it is seen.
                (*vit).do_print = ! read_seen[read_index];
                read_seen[read_index] = true;
            }
        }
    }
    else
    {
        // single end reads
        size_t current_fragment_id = SIZE_MAX;
        bool read_seen = false;

        for (vit = view; vit != view + s_last; ++vit)
        {
            if (current_fragment_id != (*vit).fragment_id)
            {
                // reset the state variables for this new fragment
                read_seen = false;
                current_fragment_id = (*vit).fragment_id;
            }
            else
            {
                // the do_print field can be initialized as soon as it is seen.
                (*vit).do_print = ! read_seen[read_index];
                (*vit).is_orphan = false;
                read_seen[read_index] = true;
            }
        }
    }

    // initialize write pointers
    char * write_pointer[2] = { fastq1_buffer_out, fastq2_buffer_out };
    char * orphan_pointer = orphan_buffer_out;

    // now, we know that all elements in this range are 'complete' and can be printed
    // to the proper buffers
    size_t nchars_printed;
    size_t read_index;
    for (vit = view; vit != view + s_last; ++vit)
    {
        if (! (*vit).do_print)
        {
            continue;
        }
        if ((*vit).is_orphan)
        {
            orphan_pointer = (*vit).print(orphan_pointer);
        }
        else
        {
            read_index = (*vit).flag.first_fragment_in_template ? 0 : 1;
            write_pointer[read_index] =  (*vit).print(write_pointer[read_index]);
        }
    }

    clock_gettime(CLOCK_REALTIME, &time_end);
    fprintf(stderr, "done. %Zu ms\n", elapsed_ms(time_begin, time_end));
    fflush(stderr);


    // clean up
    delete view;

    // the number of unprocessed lines
    return S - s_last;

}
