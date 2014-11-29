#include "sam_order.h"
#include "sam_line.h"
#include "sam_index.h"

//!!! This is muddied in that it is partially up to the user, partly a nuissance.
// SamOrder::SamOrder(SAM_ORDER _so, char const* _sort_type) : 
//     order(_so),
//     parse_fragment_id(NULL)
// {
//     omp_init_lock(&this->flowcell_hash_write_lock);

//     strcpy(this->sort_type, _sort_type);

//     switch (this->order)
//     {
//     case SAM_RID: 
//         this->less = &SamOrder::less_rid; 
//         this->equal = &SamOrder::equal_rid; 
//         break;

//     case SAM_POSITION_RID: 
//         this->less = &SamOrder::less_position_rid; 
//         this->equal = &SamOrder::equal_position_rid; 
//         break;
//     case SAM_RID_POSITION: 
//         this->less = &SamOrder::less_rid_position; 
//         this->equal = &SamOrder::equal_rid_position; 
//         break;
//     case SAM_POSITION: 
//         this->less = &SamOrder::less_position; 
//         this->equal = &SamOrder::equal_position; 
//         break;
//     case SAM_FID_POSITION: 
//         this->less = &SamOrder::less_fragment_id_position; 
//         this->equal = &SamOrder::equal_fragment_id_position; 
//         break;
//     default: 
//         fprintf(stderr, "SamOrder: unknown sort order\n");
//         exit(1);
//         break;
//     }

//     if (strcmp(this->sort_type, "MIN_ALIGN_GUIDE") == 0)
//     {
//         this->sam_index = &SamOrder::samline_position_min_align_guide;
//     }
//     else if (strcmp(this->sort_type, "ALIGN") == 0)
//     {
//         this->sam_index = &SamOrder::samline_position_align;
//     }
//     else if (strcmp(this->sort_type, "FRAGMENT") == 0)
//     {
//         this->sam_index = &SamOrder::samline_fragment_id;
//     }
//     else if (strcmp(this->sort_type, "PROJALIGN") == 0)
//     {
//         this->sam_index = &SamOrder::samline_projected_position_align;
//     }
//     else if (strcmp(this->sort_type, "NONE") == 0)
//     {
//         this->sam_index = NULL;
//     }
//     else
//     {
//         fprintf(stderr, "SamOrder: Unknown sort order: %s\n", this->sort_type);
//         exit(1);
//     }

// }


// SamOrder::SamOrder(SamOrder const& s)
// {
//     assert(false);
// }


// SamOrder::SamOrder() : order(SAM_RID), parse_fragment_id(NULL)
// {
//     sort_type[0] = '\0';
// }

// SamOrder::~SamOrder()
// {
//     // char const** keys = new char const*[contig_offsets.size()];
//     size_t k;
//     CONTIG_OFFSETS::iterator ci;
//     for (ci = contig_offsets.begin(), k = 0; ci != contig_offsets.end(); ++ci, ++k)
//     {
//         delete (*ci).first;
//         //keys[k] = (*ci).first;
//     }
//     contig_offsets.clear();

//     PROJECTIONS::const_iterator pit;
//     for (pit = this->projections.begin(); pit != this->projections.end(); ++pit)
//     {
//         delete (*pit).first;
//     }
//     this->projections.clear();

//     omp_destroy_lock(&this->flowcell_hash_write_lock);
// }


//Call to initialize the proper READ ID parsing format
// SAM_QNAME_FORMAT SamOrder::InitFromFile(FILE * sam_fh)
// {
//     SetToFirstDataLine(&sam_fh);

//     char qname[1024];

//     if (fscanf(sam_fh, "%s\t", qname) != 1)
//     {
//         qname[0] = '\0';
//     }

//     rewind(sam_fh);

//     return InitFromID(qname);

// }


// SAM_QNAME_FORMAT QNAMEFormat(char const* sam_dataline)
// {
//     SAM_QNAME_FORMAT qname_format;

//     // the settings here don't matter.
//     SamOrder dummy(SAM_RID, "ALIGN");

//     if (strlen(sam_dataline) == 0)
//     {
//         qname_format = SAM_NON_INTERPRETED;
//     }
//     else if (dummy.parse_fragment_id_numeric(sam_dataline) != QNAME_FORMAT_ERROR)
//     {
//         qname_format = SAM_NUMERIC;
//     }
//     else if (dummy.parse_fragment_id_illumina(sam_dataline) != QNAME_FORMAT_ERROR)
//     {
//         qname_format = SAM_ILLUMINA;
//     }
//     // 
//     else if (dummy.parse_fragment_id_casava_1_8(sam_dataline) != QNAME_FORMAT_ERROR)
//     {
//         qname_format = SAM_CASAVA18;
//     }
//     else
//     {
//         fprintf(stderr, "Error: SamOrder: Don't have a parser for this qname format: %s\n",
//                 sam_dataline);
//         exit(1);
//     }

//     return qname_format;
// }


// SAM_QNAME_FORMAT SamOrder::InitFromID(char const* id)
// {

//     SAM_QNAME_FORMAT qname_format = QNAMEFormat(id);

//     this->InitFromChoice(qname_format);
//     return qname_format;
// }




// void SamOrder::InitFromChoice(SAM_QNAME_FORMAT qname_format)
// {
//     switch(qname_format)
//     {
//     case SAM_NON_INTERPRETED: this->parse_fragment_id = &SamOrder::parse_fragment_id_zero; break;
//     case SAM_NUMERIC: this->parse_fragment_id = &SamOrder::parse_fragment_id_numeric; break;
//     case SAM_ILLUMINA: this->parse_fragment_id = &SamOrder::parse_fragment_id_illumina; break;
//     case SAM_CASAVA18: this->parse_fragment_id = &SamOrder::parse_fragment_id_casava_1_8; break;
//     default:
//         fprintf(stderr, "Error: InitFromChoice: Don't know this qname format\n");
//         exit(1);
//     }
// }

// bool SamOrder::Initialized() const
// {
//     return this->parse_fragment_id != NULL;
// }


//trinary comparison of two size_t's
inline int icmp(size_t a, size_t b)
{
    return a < b ? -1 : (a == b ? 0 : 1);
}



//compute the fragment id comparator
inline int fragment_id_aux(SamLine const& a, SamLine const& b)
{
    return icmp(a.fragment_id, b.fragment_id);
}


//returns the start alignment position on the catenated meta-contig,
//using contig_offsets for catenation order.
//order is consistent with [rname, pos]
size_t SamOrder::flattened_position(SamLine const* a)
{
    return samidx_flattened_position(a->rname, a->zero_based_pos(), SamOrder::contig_dictionary);
}


size_t SamOrder::flattened_position_mate(SamLine const* a) const
{
    return a->flag.is_rsam_format
        ? 0
        : samidx_flattened_position(a->next_fragment_ref_name(), a->pnext, SamOrder::contig_dictionary);
}




//for distinguishing SamLines by alignment position,
//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::less_position(SamLine const& a, SamLine const& b) const
{
    int flattened_cmp = icmp(a.flattened_pos, b.flattened_pos);

    return flattened_cmp < 0 || 
        (flattened_cmp == 0 && 
         (a.flag.this_fragment_on_neg_strand < b.flag.this_fragment_on_neg_strand ||
          (a.flag.this_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand &&
           (strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) < 0 ||
            (strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) == 0 &&
             (icmp(this->flattened_position_mate(&a),
                   this->flattened_position_mate(&b)) < 0 ||
              (icmp(this->flattened_position_mate(&a),
                    this->flattened_position_mate(&b)) == 0 &&
               a.tags.alignment_space < b.tags.alignment_space)))))));
        
}

//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::equal_position(SamLine const& a, SamLine const& b) const
{

    return 
        a.flattened_pos == b.flattened_pos &&
        a.flag.this_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand &&
        strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) == 0 &&
        flattened_position_mate(&a) == flattened_position_mate(&b) &&
        a.tags.alignment_space == b.tags.alignment_space;

}


//for distinguishing SamLines by fragment alignment position
//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::less_fposition(SamLine const& a, SamLine const& b) const
{
    size_t a_fragment_pos = a.flag.all_fragments_mapped
        ? std::min(a.flattened_pos, this->flattened_position_mate(&a))
        : a.flattened_pos;

    size_t b_fragment_pos = b.flag.all_fragments_mapped
        ? std::min(b.flattened_pos, this->flattened_position_mate(&b))
        : b.flattened_pos;

    int flattened_cmp = icmp(a_fragment_pos, b_fragment_pos);

    return flattened_cmp < 0 || 
        (flattened_cmp == 0 && 
         (a.flag.this_fragment_on_neg_strand < b.flag.this_fragment_on_neg_strand ||
          (a.flag.this_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand &&
           (strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) < 0))));

}

//consider the 'read id' to be composed of its qname and pair identity flag
//[qname, pair]
bool SamOrder::less_rid(SamLine const& a, SamLine const& b) const
{
    int qname_cmp = fragment_id_aux(a, b);
    return qname_cmp < 0 
        || (qname_cmp == 0
            && (a.flag.first_fragment_in_template < b.flag.first_fragment_in_template));
}


//[qname, pair]
bool SamOrder::equal_rid(SamLine const& a, SamLine const& b) const
{
    int qname_cmp = fragment_id_aux(a, b);
    return qname_cmp == 0
        && (a.flag.first_fragment_in_template == b.flag.first_fragment_in_template);
}


//[qname, pair, rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::less_rid_position(SamLine const& a, SamLine const& b) const
{
    return less_rid(a, b) ||
        (equal_rid(a, b) && less_position(a, b));
}

//[qname, pair, rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::equal_rid_position(SamLine const& a, SamLine const& b) const
{
    return equal_rid(a, b) && equal_position(a, b);
}


//[rname, pos, query_strand, cigar, rnext, pnext, qname, pair]
bool SamOrder::less_position_rid(SamLine const& a, SamLine const& b) const
{
    return less_position(a, b) ||
        (equal_position(a, b) && less_rid(a, b));
}

//[rname, pos, query_strand, cigar, rnext, pnext, qname, pair]
bool SamOrder::equal_position_rid(SamLine const& a, SamLine const& b) const
{
    return equal_position(a, b) && equal_rid(a, b);
}




bool SamOrder::less_fragment_id(SamLine const& a, SamLine const& b) const
{
    return fragment_id_aux(a, b) < 0;
}


bool SamOrder::equal_fragment_id(SamLine const& a, SamLine const& b) const
{
    return fragment_id_aux(a, b) == 0;
}

bool SamOrder::less_fragment_id_position(SamLine const& a, SamLine const& b) const
{
    return less_fragment_id(a, b)
        || (equal_fragment_id(a, b)
            && less_position(a, b));
}

bool SamOrder::equal_fragment_id_position(SamLine const& a, SamLine const& b) const
{
    return equal_fragment_id(a, b) && equal_position(a, b);
}






//compute the fragment id assuming a numeric read id format
// size_t SamOrder::samline_fragment_id(char const* samline)
// {
//     size_t fid = (this->*(this->parse_fragment_id))(samline);
//     if (fid == QNAME_FORMAT_ERROR)
//     {
//         fprintf(stderr, "Error: samline_fragment_id: bad qname format in this line\n%s\n\n",
//                 samline);
//         exit(1);
//     }
//     return fid;
// }

/*
less_seq_projection::less_seq_projection(SamOrder const* _so) :
    sam_order(_so) { }


bool less_seq_projection::operator()(SequenceProjection const& a,
                                     SequenceProjection const& b) const
{
    size_t a_pos = flattened_position(a.target_dna.c_str(), a.target_start_pos(),
                                      this->contig_dictionary);

    size_t b_pos = flattened_position(b.target_dna.c_str(), b.target_start_pos(), 
                                      this->contig_dictionary);

    return 
        (a_pos < b_pos
         || (a_pos == b_pos
             && (a.target_end_pos() < b.target_end_pos()
                 || (a.target_end_pos() == b.target_end_pos()
                     && (strcmp(a.source_dna.c_str(), b.source_dna.c_str()) < 0)))));
}
*/




// NON-member functions
//check that the genome contig order implied by the combination of
//transcript projections and transcript contig order is the same as
//that given by genome_contig_order
void GenerateProjectionHeader(contig_dict const* dict, FILE * out_sam_header_fh)
{
    //1. Compute a flattened genome coordinate for each transcript
    std::multimap<size_t, char const*> flat_genome_start;
    std::multimap<size_t, char const*>::const_iterator tit;

    std::map<char const*, size_t> lengths;

    // std::set<SequenceProjection>::const_iterator pit;
    // OFFSETS_ITER oit;
    HMAP::const_iterator git;
    for (size_t ord = 0; ord != dict->projected_name_map.size(); ++ord)
    {
        SequenceProjection const& sp = dict->projections[ord];

        char const* tx_contig = sp.source_dna.c_str();
        char const* tx_proj_contig = sp.target_dna.c_str();
        git = dict->name_map.find(tx_proj_contig);
        if (git == dict->name_map.end())
        {
            fprintf(stderr, "Error: GenerateProjectionHeader: didn't find genome contig %s "
                    "that was found in transcript-to-genome projection from gtf file\n",
                    tx_proj_contig);
            exit(1);
        }
        size_t genome_contig_offset = dict->offset[(*git).second];
        size_t tx_local_offset = sp.transformation.empty() ? 0 : sp.transformation[0].jump_length;
        size_t tx_length = sp.total_block_length;

    // for (pit = tx_to_genome.begin(); pit != tx_to_genome.end(); ++pit)
    // {
        // char const* tx_contig = (*pit).source_dna.c_str();
        // char const* tx_proj_contig = (*pit).target_dna.c_str();
        // oit = genome_contig_order.find(tx_proj_contig);
        // if (oit == genome_contig_order.end())
        // {
        //     fprintf(stderr, "Error: GenerateProjectionHeader: didn't find genome contig %s "
        //             "that was found in transcript-to-genome projection from gtf file\n",
        //             tx_proj_contig);
        //     exit(1);
        // }
        // size_t genome_contig_offset = (*oit).second;
        // size_t tx_local_offset = (*pit).transformation.empty() ? 0 : (*pit).transformation[0].jump_length;
        // size_t tx_length = (*pit).total_block_length;

        flat_genome_start.insert(std::make_pair(genome_contig_offset + tx_local_offset, tx_contig));
        lengths.insert(std::make_pair(tx_contig, tx_length));
    }

    fprintf(out_sam_header_fh, "@HD\tVN:samutil\tSO:coordinate\n");
    for (tit = flat_genome_start.begin(); tit != flat_genome_start.end(); ++tit)
    {
        fprintf(out_sam_header_fh, "@SQ\tSN:%s\tLN:%Zu\n", (*tit).second, lengths[(*tit).second]);
    }
}



//compute the flattened coordinate start position.  
//read id that starts out with illumina format:
//@FLOWCELL:LANE:TILE:XPOS:YPOS
// The layout, with high bits to the left, is the same
// ffffllll tttttttt ttttxxxx xxxxxxxx xxxxxxxx xxyyyyyy yyyyyyyy yyyyyyyy
// struct IlluminaID
// {
//     unsigned int ypos : 22;
//     unsigned int xpos : 22;
//     unsigned int tile : 12;
//     unsigned int lane : 4;
//     unsigned int flowcell : 4;
//     size_t get_raw() const;
// };

// size_t IlluminaID::get_raw() const
// {
//     size_t raw = static_cast<size_t>
//         (this->ypos % (1L<<22)
//          | (this->xpos % (1L<<22))<<22
//          | (this->tile % (1L<<12))<<44
//          | (this->lane % (1L<<4))<<56
//          | (this->flowcell % (1L<<4))<<60);

//     return raw;
// }


 // update record of flow cells.
 // requirements:
 /*
   Requirements:
   1. Must assign an integer in [0,16) (needs to be extended to larger range)
   2. The order in which each flowcell name is seen must NOT influence the relative order
      of any two flowcells after they are seen.
      For example, for the total set (in the proper order given) of flowcells [A,B,C,D],
      if we see them in the order [A,D,C,B], the assigned integers will be:  [8,12,10,9]
      if we see them in the order [D,C,B,A], the assigned integer will be: [8,4,2,1]
      In both cases, the final order is the same.
      However, obviously this is a problem since there are possibly many more than log2(16) = 4
      possible flowcells.  So, what to do?
   3. Once an integer as assigned to a given flowcell, it cannot be changed.
  */
// uint SamOrder::flowcell_hash_value(char const* flowcell)
// {
//     std::map<char const*, uint, less_char_ptr>::iterator fit = 
//         this->flowcell_hash.find(flowcell);

//     if (fit == this->flowcell_hash.end())
//     {
//         // need to do mutex here...
//         omp_set_lock(&this->flowcell_hash_write_lock);
//         size_t n = this->flowcell_hash.size();
//         char * fcopy = new char[strlen(flowcell) + 1];
//         strcpy(fcopy, flowcell);
//         fit = this->flowcell_hash.insert(std::make_pair(fcopy, n + 1)).first;
//         omp_unset_lock(&this->flowcell_hash_write_lock);
//     }
//     return (*fit).second;
// }
