#include "sam_helper.h"
#include "sam_to_fastq.h"
#include "time_tools.h"

#include <algorithm>
#include <parallel/algorithm>

// do not initialize anything
SAMFastqView::SAMFastqView() { }

// SAMFastqView::SAMFastqView()
//     : qname(NULL), qname_len(0), 
//       fragment_id(0), flag_raw(0), 
//       seq(NULL), qual(NULL), seqlen(0) { }

SAMFastqView::SAMFastqView(char const* _qname, size_t _qname_len, 
                           size_t _fragment_id, size_t _flag_raw, 
                           char const* _seq, char const* _qual, 
                           size_t _seqlen,
                           bool _do_print,
                           bool _is_orphan) :
    qname(_qname),
    qname_len(_qname_len),
    fragment_id(_fragment_id),
    seq(_seq),
    qual(_qual),
    seqlen(_seqlen),
    do_print(_do_print),
    is_orphan(_is_orphan)
{ 
    this->flag.set_raw(_flag_raw);
}


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
        reverse_comp_inplace(w, w + this->seqlen);
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


SAMFastqView_aux::SAMFastqView_aux(SamOrder * _sam_order) : sam_order(_sam_order) { }


// Expect a Spec-conformant samline_string, not necessarily zero-terminated
SAMFastqView SAMFastqView_aux::operator()(char * samline_string)
{
    size_t fragment_id;
    size_t flag_raw;
    char * qname = samline_string;
    char const* seq;
    char const* qual;
    int seqstart, seqend;
    int qualstart, qualend;
    size_t seqlen;
    int qname_len;

    // (taken from sam_helper.cc)
    // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
    // !!! left off here (need to find a way to
    int num_fields =
        sscanf(samline_string, 
               "%*[^\t]%n\t" "%zu\t" 
               "%*[^\t]\t" "%*u\t" "%*u\t" "%*s\t" "%*s\t" "%*u\t" "%*i\t"
               "%n%*[^\t]%n\t" "%n%*[^\t]%n\t",
               &qname_len, &flag_raw, &seqstart, &seqend, &qualstart, &qualend);

    // all but one of these fields are location fields that aren't registered in the return value count
    if (num_fields != 1)
    {
        fprintf(stderr, "Error: SAMFastqView_aux: couldn't parse samline string into a Fastq view\n");
        exit(35);
    }

    fragment_id = (this->sam_order->*(this->sam_order->sam_index))(qname);

    if (fragment_id == QNAME_FORMAT_ERROR)
    {
        fprintf(stderr, "Error: SAMFastqView_aux: invalid qname (not Illumina format): %s\n", qname);
    }
    seq = samline_string + seqstart;
    qual = samline_string + qualstart;
    seqlen = seqend - seqstart;

    return SAMFastqView(qname, qname_len, fragment_id, flag_raw, seq, qual, seqlen, false, false);
}


void init_fastq_view(std::vector<char *> & sam_lines,
                     SAMFastqView * fastq_view)
{
    if (sam_lines.empty())
    {
        return;
    }
    // We are simply taking advantage of sam_order's facility to parse the qname field
    // and transform it into a fragment_id index.  The "FRAGMENT" part is of no consequence.

    SamOrder sam_order(SAM_RID, "FRAGMENT");
    sam_order.InitFromID(sam_lines[0]);

    SAMFastqView_aux s2v_aux(& sam_order);
    __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), fastq_view, s2v_aux);

}


// finds the upper bound in fastq_view which is the highest record
// that is preceded by a record of a different fragment type, i.e.
// in [A, A, A, B, B, C, C, C, {D}, D], finds the bracketed record {D}
// and returns its index
size_t find_last_fragment_bound(SAMFastqView * fastq_view, size_t S, bool is_last_chunk)
{
    if (S == 0 || is_last_chunk)
    {
        return S;
    }
    size_t s_last = S - 1;
    size_t last_fragment_id = fastq_view[s_last].fragment_id;
    
    while(s_last > 0 && fastq_view[s_last].fragment_id == last_fragment_id)
    {
        --s_last;
    }
    return s_last + 1;
}                     

// post-condition: any non-null member of marked_reads
// will have its is_orphan field correctly set
void init_orphan_field(SAMFastqView * marked_reads[])
{
    bool is_orphan = (marked_reads[0] == NULL || marked_reads[1] == NULL);
    // assert(is_orphan == false);
    // the is_orphan field must wait until seeing all reads of this fragment
    if (marked_reads[0] != NULL)
    {
        (*marked_reads[0]).is_orphan = is_orphan;
    }
    if (marked_reads[1] != NULL)
    {
        (*marked_reads[1]).is_orphan = is_orphan;
    }
}


// build a fastq 'view' of the SAM records, starting from
// an initialized fastq_view (initialized with init_fastq_view.
// this view stores pointers to the relevant sections of the SAM line
// and also flags for do_print and is_orphan
void set_flags_fastq_view_paired(std::vector<char *> & sam_lines,
                                 SAMFastqView * fastq_view,
                                 bool is_last_chunk)
{
    size_t S = sam_lines.size();

    // 2. traverse the index in a serial loop, dealing and revcomp'ing.
    // find the bound at which we can guarantee complete entries
    size_t s_last = find_last_fragment_bound(fastq_view, S, is_last_chunk);
    SAMFastqView *vit_end = fastq_view + s_last;
        
    SAMFastqView * marked_reads[2] = { NULL, NULL };
    size_t read_index; // 0 or 1 depending on whether first in template or not.

    // current_fragment_id is initialized so that the first iteration of the loop
    // looks like it is NOT a new chunk
    size_t current_fragment_id = (s_last > 0) ? fastq_view[0].fragment_id : SIZE_MAX;
    for (SAMFastqView * vit = fastq_view; vit != vit_end; ++vit)
    {
        read_index = (*vit).flag.first_fragment_in_template ? 0 : 1;

        // update the is_orphan state of the marked reads by looking at the pair.
        // if we are in a new fragment, initialize state variables to their naive state.
        if (current_fragment_id != (*vit).fragment_id)
        {
            init_orphan_field(marked_reads);

            // reset the state variables for this new fragment
            marked_reads[0] = NULL;
            marked_reads[1] = NULL;
        }

        // the do_print field can be initialized as soon as it is seen.
        if (marked_reads[read_index] == NULL)
        {
            // we haven't seen a read of this index -- this is the one we will print
            (*vit).do_print = true;
            marked_reads[read_index] = vit;
        }
        current_fragment_id = (*vit).fragment_id;
    }

    // this is necessary at the end of the loop
    init_orphan_field(marked_reads);
   
}


// initialize the do_print flag of an initialized fastq_view
void set_flags_fastq_view_single(std::vector<char *> & sam_lines,
                                 SAMFastqView * fastq_view,
                                 bool is_last_chunk)
{
    size_t S = sam_lines.size();
    size_t s_last = find_last_fragment_bound(fastq_view, S, is_last_chunk);
    SAMFastqView *vit_end = fastq_view + s_last;

    // single end reads
    SAMFastqView * marked_read = NULL;

    size_t current_fragment_id = (s_last > 0) ? fastq_view[0].fragment_id : SIZE_MAX;
    for (SAMFastqView * vit = fastq_view; vit != vit_end; ++vit)
    {
        if (current_fragment_id != (*vit).fragment_id)
        {
            // reset the state variables for this new fragment
            marked_read = NULL;
        }

        // the do_print field can be initialized as soon as it is seen.
        (*vit).do_print = (marked_read == NULL);
        marked_read = vit;
        current_fragment_id = (*vit).fragment_id;
    }
}


// writes vS fastq-formatted paired-end reads to the appropriate
// buffers in uncompressed format.  Records the number of bytes written
// to each of the three buffers in the output parameter nbytes_written
void write_fastq_view_paired(SAMFastqView * fastq_view,
                             size_t vS,
                             char * fastq1_buffer_out,
                             char * fastq2_buffer_out,
                             char * orphan_buffer_out,
                             size_t * nbytes_written)
{
    
    // initialize write pointers
    char * write_pointer[2] = { fastq1_buffer_out, fastq2_buffer_out };
    char * orphan_pointer = orphan_buffer_out;
    
    // now, we know that all elements in this range are 'complete' and
    // can be printed to the proper buffers
    size_t read_index;
    SAMFastqView * vit_end = fastq_view + vS;
    for (SAMFastqView * vit = fastq_view; vit != vit_end; ++vit)
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
            write_pointer[read_index] = (*vit).print(write_pointer[read_index]);
        }
    }
    nbytes_written[0] = (write_pointer[0] - fastq1_buffer_out);
    nbytes_written[1] = (write_pointer[1] - fastq2_buffer_out);
    nbytes_written[2] = (orphan_pointer - orphan_buffer_out);
    // printf("nbytes_written[2] = %Zu\n", nbytes_written[2]);
    // printf("orphan_pointer - orphan_buffer_out = %li\n", (orphan_pointer - orphan_buffer_out));
}


// writes vS fastq-formatted single-end reads to the fastq_buffer_out
// in uncompressed format.  Returns the number of bytes written
void write_fastq_view_single(SAMFastqView * fastq_view,
                             size_t vS,
                             char * fastq_buffer_out,
                             size_t * nbytes_written)
{
    // initialize write pointers
    char * write_pointer = fastq_buffer_out;
    
    SAMFastqView * vit_end = fastq_view + vS;
    for (SAMFastqView * vit = fastq_view; vit != vit_end; ++vit)
    {
        if ((*vit).do_print)
        {
            write_pointer = (*vit).print(write_pointer);
        }
    }
    *nbytes_written = (write_pointer - fastq_buffer_out);
}
