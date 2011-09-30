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

SEQ_PROJ_ITER 
non_overlapping_range(SEQ_PROJ_ITER start, SEQ_PROJ_ITER set_end)
{
    SEQ_PROJ_ITER end = start;
    SEQ_PROJ_ITER max = start;
    
    if (start == set_end)
    {
        return start;
    }

    while (strcmp((*start).target_dna.c_str(), (*end).target_dna.c_str()) == 0
           && ((*end).target_start_pos() < (*start).target_end_pos())
           && end != set_end)
    {
        max = ((*end).target_end_pos() > (*max).target_end_pos())
            ? end
            : max;
        ++end;
    }
    
    if (max != start)
    {
        return non_overlapping_range(max, set_end);
    }
    else
    {
        return end;
    }
}




/*
  returns true if the named contigs do not have overlapping
  projections on the genome. by definition, unmapped contigs have no
  projection, and are thus non-overlapping.
 */
/*
bool tx_nonoverlapping_or_passthrough(PROJ_MAP const& tx_projections,
                                      char const* tx1, char const* tx2)
{
    if (strcmp(tx1, "*") == 0 || strcmp(tx2, "*") == 0)
    {
        return true;
    }
    else if (strcmp(tx1, tx2) == 0)
    {
        return false;
    }

    // check if these are non-overlapping contigs
    PROJ_MAP::const_iterator proj1_iter = tx_projections.find(tx1);
    PROJ_MAP::const_iterator proj2_iter = tx_projections.find(tx2);

    if (proj1_iter == tx_projections.end() || proj2_iter == tx_projections.end())
    {
        //either of these are not projected, so we're in pass-through mode.
        return true;
    }
                    
    SequenceProjection const& sp1 = *(*proj1_iter).second;
    SequenceProjection const& sp2 = *(*proj2_iter).second;

    return 
        sp1.target_dna != sp2.target_dna
        || sp1.target_end_pos() <= sp2.target_start_pos()
        || sp2.target_end_pos() <= sp1.target_start_pos();
}
*/


project_dedup_print::project_dedup_print(PROJ_MAP const* _proj_map,
                                         SamOrder const* _tx_ord,
                                         SamOrder const* _genome_ord,
                                         bool _inserts_are_introns,
                                         bool _retain_unseq,
                                         size_t _avg_len)
    : projections(_proj_map), 
      tx_sam_order(_tx_ord),
      genome_sam_order(_genome_ord),
      inserts_are_introns(_inserts_are_introns),
      retain_unsequenced_projection(_retain_unseq),
      average_rsam_line_length(_avg_len) { }


std::vector<char> * 
project_dedup_print::operator()(std::pair<SAMIT, SAMIT> const& range)
{
    SamBuffer output_buf(this->genome_sam_order, false);

    // project them
    SINGLE_READ_SET::iterator beg, sit, end;

    PROJ_MAP::const_iterator proj_iter;

    char last_contig[256] = "";
    CONTIG_OFFSETS::const_iterator contig_iter = 
        this->genome_sam_order->contig_offsets.begin();

    char transcript_strand[2] = "";

    SAMIT cur;
    for (cur = range.first; cur != range.second; ++cur)
    {
        SamLine * samline = const_cast<SamLine *>(*cur);

        // here is where we should efficiently compare and replace / ignore based on
        // mapping quality.
        if (strcmp(last_contig, (*cur)->rname) != 0)
        {
            // update proj_iter if necessary
            proj_iter = this->projections->find((*cur)->rname);
            transcript_strand[0] = proj_iter == this->projections->end()
                ? '.' 
                : ((*(*proj_iter).second).same_strand ? '+' : '-');
        }

        if ((*cur)->flag.all_fragments_mapped)
        {
            if (proj_iter != this->projections->end())
            {
                ApplySequenceProjection(*(*proj_iter).second, samline, 
                                        this->inserts_are_introns);

                samline->add_tag("XS", 'A', transcript_strand);
            }
        }
        
        samline->SetFlattenedPosition(this->genome_sam_order->contig_offsets, 
                                      & contig_iter);

        if (! this->retain_unsequenced_projection)
        {
            samline->SetCigarCompared();
        }

        InsertResult ins = output_buf.insert(samline);
        
        strcpy(last_contig, (*ins.surviving_entry)->rname);

        if (! ins.was_inserted)
        {
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


    }

    size_t num_entries = output_buf.unique_entries.size();

    std::vector<char> * out_lines = new std::vector<char>;
    out_lines->reserve(num_entries * this->average_rsam_line_length);
    
    char tmp_buffer[1024];
    //char * write_pointer = out_buffer;
    for (sit = output_buf.unique_entries.begin(); 
         sit != output_buf.unique_entries.end(); ++sit)
    {
        (*sit)->sprint(tmp_buffer);
        delete (*sit);
        out_lines->insert(out_lines->end(), tmp_buffer, tmp_buffer + strlen(tmp_buffer));
    }
    return out_lines;
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
        alloc_buf = new char[strlen(tmp_buf) + 1];
        strcpy(alloc_buf, tmp_buf);
    }
    delete samline;

    return alloc_buf;
}
