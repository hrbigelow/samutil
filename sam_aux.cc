#include "sam_aux.h"


/*
  pseudo-code:

  precondition: start_iter is a valid start bound.  (i.e. it is at the
  beginning of tx_projections, or the previous one does not overlap in
  genome extent.

  precondition: tx_projections are sorted by contig, contig_start.

  end_iter = start_iter;
  iter max_iter = end_iter;
  while (end_iter.contig == start_iter.contig 
         && end_iter.start < start_iter.end
         && end_iter != tx_projections.end())
  {
       max_iter = end_iter.end > max_iter.end ? end_iter : max_iter;
       ++end_iter;
  }
  if (max_iter != start_iter)
  {
     // recurse on max_iter
  }
  else
  {
     // return end_iter
  }    

  postcondition: returned value is a valid bound, such that [start, end) and [end, ...)
  do not overlap in extent on genome.

*/


// Find a breakpoint where  
/*
SEQ_PROJ_ITER 
non_overlapping_range(SEQ_PROJ_ITER start, SEQ_PROJ_ITER set_end)
{
    SEQ_PROJ_ITER end = start;
    SEQ_PROJ_ITER max = start;
    
    if (start == set_end)
    {
        return start;
    }

    // find the projection 'max' that overlaps 'start' and has the
    // furthest end_pos(). this projection may indeed *be* 'start'
    while (strcmp((*start).target_dna.c_str(), (*end).target_dna.c_str()) == 0
           && ((*end).target_start_pos() < (*start).target_end_pos())
           && end != set_end)
    {
        max = ((*max).target_end_pos() < (*end).target_end_pos())
            ? end
            : max;
        ++end;
    }
    
    if (max != start)
    {
        // recurse, since now there may be another max that overlaps
        // this max.
        return non_overlapping_range(max, set_end);
    }
    else
    {
        return end;
    }
}
*/





pdp_result::pdp_result() :
    lines(NULL),
    num_records_retained(0),
    num_records_discarded(0) { }


project_dedup_print::project_dedup_print(SamOrder const* _genome_ord,
                                         bool _inserts_are_introns,
                                         bool _retain_unseq,
                                         size_t _avg_len)
    : genome_sam_order(_genome_ord),
      inserts_are_introns(_inserts_are_introns),
      retain_unsequenced_projection(_retain_unseq),
      average_rsam_line_length(_avg_len) { }


pdp_result
project_dedup_print::operator()(std::pair<SAMIT, SAMIT> const& range)
{
    SamBuffer output_buf(this->genome_sam_order, false);
    contig_dict * dict = this->contig_dictionary;

    pdp_result result;

    // project them
    SINGLE_READ_SET::iterator beg, sit, end;

    // PROJECTIONS::const_iterator proj_iter = this->genome_sam_order->projections.end();
    HMAP::const_iterator proj_iter = dict->projected_name_map.end();

    char prev_unprojected_rname[256] = "";

    // CONTIG_OFFSETS::const_iterator contig_iter = 
    //     this->genome_sam_order->contig_offsets.begin();
    HMAP::const_iterator contig_iter = dict->name_map.begin();

    char transcript_strand[2] = "";

    SAMIT cur;
    for (cur = range.first; cur != range.second; ++cur)
    {
        SamLine * samline = const_cast<SamLine *>(*cur);

        // here is where we should efficiently compare and replace /
        // ignore based on mapping quality.
        if (strcmp(prev_unprojected_rname, (*cur)->rname) != 0)
        {
            // update proj_iter if necessary
            // proj_iter = this->genome_sam_order->projections.find((*cur)->rname);
            // transcript_strand[0] = proj_iter == this->genome_sam_order->projections.end()
            //     ? '.' 
            //     : ((*proj_iter).second.same_strand ? '+' : '-');

            proj_iter = dict->projected_name_map.find((*cur)->rname);
            transcript_strand[0] = (proj_iter == dict->projected_name_map.end())
                ? '.'
                : (dict->projections[(*proj_iter).second].same_strand ? '+' : '-');

            strcpy(prev_unprojected_rname, (*cur)->rname);
        }

        if ((*cur)->flag.all_fragments_mapped)
        {
            // if (proj_iter != this->genome_sam_order->projections.end())
            // {
            //     ApplySequenceProjection((*proj_iter).second, samline, this->inserts_are_introns);
            //     samline->add_tag("XS", 'A', transcript_strand);
            // }
            if (proj_iter != dict->projected_name_map.end())
            {
                ApplySequenceProjection(dict->projections[(*proj_iter).second],
                                        samline,
                                        this->inserts_are_introns);
                samline->add_tag("XS", 'A', transcript_strand);
            }

        }
        
        // samline->SetFlattenedPosition(this->genome_sam_order->contig_offsets, 
        // & contig_iter);
        samline->SetFlattenedPosition(contig_dictionary);

        if (! this->retain_unsequenced_projection)
        {
            samline->SetCigarCompared();
        }

        InsertResult ins = output_buf.insert(samline);
        
        if (! ins.was_inserted)
        {
            ++result.num_records_discarded;

            // there is a collision.  now we must prioritize.
            if (ins.remaining_entry->mapq > (*ins.surviving_entry)->mapq)
            {
                output_buf.replace(ins.surviving_entry, ins.remaining_entry);
            }
            else
            {
                delete ins.remaining_entry;
            }
        }
        else
        {
            ++result.num_records_retained;
        }

    }

    size_t num_entries = output_buf.unique_entries.size();

    result.lines = new std::vector<char>;
    result.lines->reserve(num_entries * this->average_rsam_line_length);
    
    char tmp_buffer[1024];
    //char * write_pointer = out_buffer;
    for (sit = output_buf.unique_entries.begin(); 
         sit != output_buf.unique_entries.end(); ++sit)
    {
        (*sit)->sprint(tmp_buffer);
        delete (*sit);
        result.lines->insert(result.lines->end(), tmp_buffer, tmp_buffer + strlen(tmp_buffer));
    }
    return result;
}


rsam_to_sam_binary::rsam_to_sam_binary(char const* _seq_buffer,
                                       size_t _data_buffer_offset,
                                       SamFilter const* _sf)
    : seq_buffer(_seq_buffer), 
      data_buffer_offset(_data_buffer_offset),
      sam_filter(_sf)
{ }
                                       

// convert an rsam-formatted line into a string holding the printed
// set of SAM records, with joined QNAME, SEQ, and QUAL fields
char * rsam_to_sam_binary::operator()(SamLine * samline, INDEX_ITER li_iter)
{
    
    assert((*li_iter).start_offset >= this->data_buffer_offset);
    char * alloc_buf = NULL;

    if (this->sam_filter->pass(samline))
    {
        // find the zero-terminated seq_data from this iter and an offset
        char const* seq_data = 
            this->seq_buffer + (*li_iter).start_offset - this->data_buffer_offset;

        char tmp_buf[4096 * 8];
        samline->print_rsam_as_sam(seq_data, tmp_buf);
        assert(strlen(tmp_buf) < 4096 * 7);

        alloc_buf = new char[strlen(tmp_buf) + 1];
        strcpy(alloc_buf, tmp_buf);
    }

    return alloc_buf;
}


unmapped_rsam_to_fastq::unmapped_rsam_to_fastq(char const* _seq_buffer,
                                               size_t _data_buffer_offset)
    : seq_buffer(_seq_buffer), 
      data_buffer_offset(_data_buffer_offset)
{ }
                                       

// convert an rsam-formatted line into a string holding the printed
// set of SAM records, with joined QNAME, SEQ, and QUAL fields
char const* unmapped_rsam_to_fastq::operator()(SamLine * samline, INDEX_ITER li_iter)
{
    
    assert((*li_iter).start_offset >= this->data_buffer_offset);
    char const* seq_ptr = NULL;

    if (! samline->flag.all_fragments_mapped
        && samline->tags.stratum_rank == 1)
    {
        // find the zero-terminated seq_data from this iter and an offset
        seq_ptr = this->seq_buffer + (*li_iter).start_offset - this->data_buffer_offset;
        assert(strlen(seq_ptr) > 0);
    }

    return seq_ptr;
}


delete_samline::delete_samline() { }

void delete_samline::operator()(SamLine * samline)
{
    delete samline;
}


// assume a tab-separated line of [id seq1 qual1 [seq2 qual2 [...]]]
// line is null-terminated
void print_fqd_as_fastq(char const* fqd_line, FILE * fq_fh)
{
    char const* seq_start = strchr(fqd_line, '\t');
    size_t id_length = std::distance(fqd_line, seq_start);

    while (*seq_start++ != '\n')
    {
        fputs("@", fq_fh);
        fwrite(fqd_line, 1, id_length, fq_fh);
        fputs("\n", fq_fh);
        // seq string
        seq_start += fwrite(seq_start, 1, strcspn(seq_start, "\t"), fq_fh) + 1;
        fputs("\n+\n", fq_fh);
        // qual string; will increment seq_start past newline
        seq_start += fwrite(seq_start, 1, strcspn(seq_start, "\t\n"), fq_fh);
        fputs("\n", fq_fh);
    }
}


// find the average length of the printed string representation of
// SamLines
size_t average_line_length(std::vector<SamLine *>::const_iterator beg,
                           std::vector<SamLine *>::const_iterator end,
                           size_t nsample)
{
    size_t average_length = 200; // reasonable value in absence of any other information
    size_t sum_of_lengths = 0;

    std::vector<SamLine *>::const_iterator it_end = beg;
    std::advance(it_end, std::min(static_cast<size_t>(std::distance(beg, end)), nsample));

    char write_buffer[1024];
    for ( ; beg != it_end; ++beg)
    {
        (*beg)->sprint(write_buffer);
        sum_of_lengths += strlen(write_buffer);
    }
    if (nsample > 0)
    {
        average_length = sum_of_lengths / nsample;
    }
    return average_length;
}
