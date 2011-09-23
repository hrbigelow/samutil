#include "sam_aux.h"

/*
  returns true if the named contigs do not have overlapping
  projections on the genome. by definition, unmapped contigs have no
  projection, and are thus non-overlapping.
 */
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


project_dedup_print::project_dedup_print(PROJ_MAP const* _proj_map,
                                         SamOrder const* _tx_ord,
                                         SamOrder const* _genome_ord,
                                         bool _inserts_are_introns,
                                         size_t _avg_len)
    : projections(_proj_map), 
      tx_sam_order(_tx_ord),
      genome_sam_order(_genome_ord),
      inserts_are_introns(_inserts_are_introns),
      average_rsam_line_length(_avg_len) { }


std::vector<char> * 
project_dedup_print::operator()(std::pair<SAMIT, SAMIT> const& range)
{
    // SamBuffer input_buf(this->tx_sam_order, false);
    SamBuffer output_buf(this->genome_sam_order, false);

    // SAMIT cur;
    // for (cur = range.first; cur != range.second; ++cur)
    // {
    //     input_buf.insert(*cur);
    // }
    // // At this point, the buffer should have joined everything together
    // if (! input_buf.incomplete_entries.empty())
    // {
    //     fprintf(stderr, "Error: at a fragment boundary but there are still incomplete SAM entries\n");
    //     exit(1);
    // }

    // project them
    SINGLE_READ_SET::iterator beg, sit, end;

    PROJ_MAP::const_iterator proj_iter;

    char last_contig[256] = "";
    CONTIG_OFFSETS::const_iterator contig_iter;

    // for (sit = input_buf.unique_entries.begin(); 
    //      sit != input_buf.unique_entries.end(); ++sit)
    SAMIT cur;
    for (cur = range.first; cur != range.second; ++cur)
    {
        SamLine * samline = const_cast<SamLine *>(*cur);

        if (strcmp(last_contig, (*cur)->rname) != 0)
        {
            contig_iter = this->tx_sam_order->contig_offsets.find((*cur)->rname);
            assert(contig_iter != this->tx_sam_order->contig_offsets.end());

            if ((*cur)->flag.all_fragments_mapped)
            {
                proj_iter = this->projections->find((*cur)->rname);
                if (proj_iter != this->projections->end())
                {
                    ApplySequenceProjection(*(*proj_iter).second, samline, 
                                            this->inserts_are_introns);
                }
            }
            strcpy(last_contig, (*cur)->rname);

        }

        if (output_buf.insert(samline).second)
        {
            samline->SetFlattenedPosition(this->genome_sam_order->contig_offsets, & contig_iter);
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
                                       size_t _data_buffer_offset)
    : seq_buffer(_seq_buffer), 
      data_buffer_offset(_data_buffer_offset) { }
                                       

// convert an rsam-formatted line into a string holding the printed
// set of SAM records, with joined QNAME, SEQ, and QUAL fields
char * rsam_to_sam_binary::operator()(SamLine * samline, INDEX_ITER li_iter)
{
    char tmp_buf[4096 * 8];
    
    assert((*li_iter).start_offset >= this->data_buffer_offset);

    // find the zero-terminated seq_data from this iter and an offset
    char const* seq_data = 
        this->seq_buffer + (*li_iter).start_offset - this->data_buffer_offset;

    samline->print_rsam_as_sam(seq_data, tmp_buf);

    char * alloc_buf = new char[strlen(tmp_buf) + 1];
    strcpy(alloc_buf, tmp_buf);

    delete samline;

    return alloc_buf;
}
