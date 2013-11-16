#include "sam_order.h"
#include "sam_helper.h"
#include "file_utils.h"
#include "gtf.h"

#include <set>
#include <omp.h>

//!!! This is muddied in that it is partially up to the user, partly a nuissance.
SamOrder::SamOrder(SAM_ORDER _so, char const* _sort_type) : 
    order(_so),
    parse_fragment_id(NULL)
{
    omp_init_lock(&this->flowcell_hash_write_lock);

    strcpy(this->sort_type, _sort_type);

    switch (this->order)
    {
    case SAM_RID: 
        this->less = &SamOrder::less_rid; 
        this->equal = &SamOrder::equal_rid; 
        break;

    case SAM_POSITION_RID: 
        this->less = &SamOrder::less_position_rid; 
        this->equal = &SamOrder::equal_position_rid; 
        break;
    case SAM_RID_POSITION: 
        this->less = &SamOrder::less_rid_position; 
        this->equal = &SamOrder::equal_rid_position; 
        break;
    case SAM_POSITION: 
        this->less = &SamOrder::less_position; 
        this->equal = &SamOrder::equal_position; 
        break;
    case SAM_FID_POSITION: 
        this->less = &SamOrder::less_fragment_id_position; 
        this->equal = &SamOrder::equal_fragment_id_position; 
        break;
    default: 
        fprintf(stderr, "SamOrder: unknown sort order\n");
        exit(1);
        break;
    }

    if (strcmp(this->sort_type, "MIN_ALIGN_GUIDE") == 0)
    {
        this->sam_index = &SamOrder::samline_position_min_align_guide;
    }
    else if (strcmp(this->sort_type, "ALIGN") == 0)
    {
        this->sam_index = &SamOrder::samline_position_align;
    }
    else if (strcmp(this->sort_type, "FRAGMENT") == 0)
    {
        this->sam_index = &SamOrder::samline_fragment_id;
    }
    else if (strcmp(this->sort_type, "PROJALIGN") == 0)
    {
        this->sam_index = &SamOrder::samline_projected_position_align;
    }
    else if (strcmp(this->sort_type, "NONE") == 0)
    {
        this->sam_index = NULL;
    }
    else
    {
        fprintf(stderr, "SamOrder: Unknown sort order: %s\n", this->sort_type);
        exit(1);
    }

}


// SamOrder::SamOrder(SamOrder const& s)
// {
//     assert(false);
// }


SamOrder::SamOrder() : order(SAM_RID), parse_fragment_id(NULL)
{
    sort_type[0] = '\0';
}

SamOrder::~SamOrder()
{
    // char const** keys = new char const*[contig_offsets.size()];
    size_t k;
    CONTIG_OFFSETS::iterator ci;
    for (ci = contig_offsets.begin(), k = 0; ci != contig_offsets.end(); ++ci, ++k)
    {
        delete (*ci).first;
        //keys[k] = (*ci).first;
    }
    contig_offsets.clear();

    PROJECTIONS::const_iterator pit;
    for (pit = this->projections.begin(); pit != this->projections.end(); ++pit)
    {
        delete (*pit).first;
    }
    this->projections.clear();

    omp_destroy_lock(&this->flowcell_hash_write_lock);
}


//Call to initialize the proper READ ID parsing format
SAM_QNAME_FORMAT SamOrder::InitFromFile(FILE * sam_fh)
{
    SetToFirstDataLine(&sam_fh);

    char qname[1024];

    if (fscanf(sam_fh, "%s\t", qname) != 1)
    {
        qname[0] = '\0';
    }

    rewind(sam_fh);

    return InitFromID(qname);

}


SAM_QNAME_FORMAT QNAMEFormat(char const* sam_dataline)
{
    SAM_QNAME_FORMAT qname_format;

    // the settings here don't matter.
    SamOrder dummy(SAM_RID, "ALIGN");

    if (strlen(sam_dataline) == 0)
    {
        qname_format = SAM_NON_INTERPRETED;
    }
    else if (dummy.parse_fragment_id_numeric(sam_dataline) != QNAME_FORMAT_ERROR)
    {
        qname_format = SAM_NUMERIC;
    }
    else if (dummy.parse_fragment_id_illumina(sam_dataline) != QNAME_FORMAT_ERROR)
    {
        qname_format = SAM_ILLUMINA;
    }
    // 
    else if (dummy.parse_fragment_id_casava_1_8(sam_dataline) != QNAME_FORMAT_ERROR)
    {
        qname_format = SAM_CASAVA18;
    }
    else
    {
        fprintf(stderr, "Error: SamOrder: Don't have a parser for this qname format: %s\n",
                sam_dataline);
        exit(1);
    }

    return qname_format;
}


SAM_QNAME_FORMAT SamOrder::InitFromID(char const* id)
{

    SAM_QNAME_FORMAT qname_format = QNAMEFormat(id);

    this->InitFromChoice(qname_format);
    return qname_format;
}




void SamOrder::InitFromChoice(SAM_QNAME_FORMAT qname_format)
{
    switch(qname_format)
    {
    case SAM_NON_INTERPRETED: this->parse_fragment_id = &SamOrder::parse_fragment_id_zero; break;
    case SAM_NUMERIC: this->parse_fragment_id = &SamOrder::parse_fragment_id_numeric; break;
    case SAM_ILLUMINA: this->parse_fragment_id = &SamOrder::parse_fragment_id_illumina; break;
    case SAM_CASAVA18: this->parse_fragment_id = &SamOrder::parse_fragment_id_casava_1_8; break;
    default:
        fprintf(stderr, "Error: InitFromChoice: Don't know this qname format\n");
        exit(1);
    }
}


void SamOrder::InitProjection(char const* gtf_file)
{
    char const* species = "dummy";
    std::set<SequenceProjection> tx_to_genome = 
        gtf_to_sequence_projection(gtf_file, species);

    std::set<SequenceProjection>::const_iterator tgi;
    for (tgi = tx_to_genome.begin(); tgi != tx_to_genome.end(); ++tgi)
    {
        char * transcript = new char[(*tgi).source_dna.size() + 1];
        strcpy(transcript, (*tgi).source_dna.c_str());
        this->projections.insert(std::make_pair(transcript, SequenceProjection(*tgi)));
    }
}

bool SamOrder::Initialized() const
{
    return this->parse_fragment_id != NULL;
}


//trinary comparison of two size_t's
inline int icmp(size_t a, size_t b)
{
    return a < b ? -1 : (a == b ? 0 : 1);
}

size_t flattened_position_aux(char const* contig,
                              size_t position, 
                              CONTIG_OFFSETS const& contig_offsets,
                              CONTIG_OFFSETS::const_iterator * contig_iter)
{
    
    size_t index;
    *contig_iter = contig_offsets.find(contig);
        
    if (*contig_iter == contig_offsets.end())
    {
        fprintf(stderr, "flattened_position_aux: error: rname %s (at %Zu) "
                "does not exist in provided contig index\n",
                contig, position);
        exit(1);
    }
    
    index = (**contig_iter).second + position;
    
    return index;
}


//compute the fragment id comparator
inline int fragment_id_aux(SamLine const& a, SamLine const& b)
{
    return icmp(a.fragment_id, b.fragment_id);
}


//returns the start alignment position on the catenated meta-contig,
//using contig_offsets for catenation order.
//order is consistent with [rname, pos]
size_t SamOrder::flattened_position(SamLine const* a, 
                                    CONTIG_OFFSETS::const_iterator * contig_iter) const
{
    
    return flattened_position_aux(a->rname, a->zero_based_pos(),
                                  this->contig_offsets,
                                  contig_iter);
}


size_t SamOrder::flattened_position_mate(SamLine const* a, 
                                         CONTIG_OFFSETS::const_iterator * contig_iter) const
{
    return a->flag.is_rsam_format
        ? 0
        : flattened_position_aux(a->next_fragment_ref_name(), a->pnext,
                                 this->contig_offsets,
                                 contig_iter);
}




//for distinguishing SamLines by alignment position,
//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::less_position(SamLine const& a, SamLine const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;
    int flattened_cmp = icmp(a.flattened_pos, b.flattened_pos);
    // int flattened_cmp = icmp(this->flattened_position(&a, &dummy),
    //                          this->flattened_position(&b, &dummy));

    return flattened_cmp < 0 || 
        (flattened_cmp == 0 && 
         (a.flag.this_fragment_on_neg_strand < b.flag.this_fragment_on_neg_strand ||
          (a.flag.this_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand &&
           (strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) < 0 ||
            (strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) == 0 &&
             (icmp(this->flattened_position_mate(&a, &dummy),
                   this->flattened_position_mate(&b, &dummy)) < 0 ||
              (icmp(this->flattened_position_mate(&a, &dummy),
                    this->flattened_position_mate(&b, &dummy)) == 0 &&
               a.tags.alignment_space < b.tags.alignment_space)))))));
        
}

//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::equal_position(SamLine const& a, SamLine const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;

    return 
        a.flattened_pos == b.flattened_pos &&
        a.flag.this_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand &&
        strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) == 0 &&
        flattened_position_mate(&a, &dummy) == flattened_position_mate(&b, &dummy) &&
        a.tags.alignment_space == b.tags.alignment_space;

    // flattened_position(&a, &dummy) == flattened_position(&b, &dummy) &&
    // a.this_fragment_on_pos_strand() == b.this_fragment_on_pos_strand() &&
    // strcmp(a.cigar_for_comparison(), b.cigar_for_comparison()) == 0 &&
    // flattened_position_mate(&a, &dummy) == flattened_position_mate(&b, &dummy);
}


//for distinguishing SamLines by fragment alignment position
//[rname, pos, query_strand, cigar, rnext, pnext]
bool SamOrder::less_fposition(SamLine const& a, SamLine const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;

    size_t a_fragment_pos = a.flag.all_fragments_mapped
        ? std::min(a.flattened_pos, this->flattened_position_mate(&a, &dummy))
        : a.flattened_pos;

    size_t b_fragment_pos = b.flag.all_fragments_mapped
        ? std::min(b.flattened_pos, this->flattened_position_mate(&b, &dummy))
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




//compute contig lengths from sam_fh.
//rewinds sam_fh to beginning
void SamOrder::AddHeaderContigStats(char const* header)
{
    size_t contig_offset = 0;
    char * hdr_copy = new char[strlen(header) + 1];
    strcpy(hdr_copy, header);

    char * last_fragment;
    std::vector<char *> header_lines =
        FileUtils::find_complete_lines_nullify(hdr_copy, & last_fragment);
    
    assert(strlen(last_fragment) == 0);

    char const* line;
    //while (! feof(sam_fh))
    for (std::vector<char *>::const_iterator it = header_lines.begin();
         it != header_lines.end(); ++it)
    {
        //if there is a SQ tag, parse and record
        line = *it;
        if (strncmp(line, "@SQ", 3) == 0)
        {
            //this is a SQ line
            char contig_name[1000];
            
            size_t contig_length;
            // need to convert this to something tag-agnostic
            sscanf(line, "@SQ\tSN:%[^\t]\tLN:%zu", contig_name, &contig_length);
            this->contig_lengths[std::string(contig_name)] = contig_length;

            char * contig_name_copy = new char[strlen(contig_name) + 1];
            strcpy(contig_name_copy, contig_name);

            this->contig_offsets[contig_name_copy] = contig_offset;
            contig_offset += contig_length;
        }
    }
    delete hdr_copy;

    this->contig_lengths[std::string("*")] = 0;

    
    char * contig_name_copy = new char[2];
    strcpy(contig_name_copy, "*");

    this->contig_offsets[contig_name_copy] = contig_offset;

    if (contig_offset == 0)
    {
        fprintf(stderr, "No SQ lines in header, or no SN / LN fields defining contigs\n");
        exit(1);
    }

    this->contig_offsets.rehash(this->contig_offsets.size() * 10);

    // fprintf(stderr, "size()=%zu, bucket_count()=%zu, load_factor()=%f\n", 
    //         this->contig_offsets.size(),
    //         this->contig_offsets.bucket_count(),
    //         this->contig_offsets.load_factor());

    std::map<size_t, size_t> bucket_hist;
    for (size_t i = 0; i != this->contig_offsets.bucket_count(); ++i)
    {
        bucket_hist[this->contig_offsets.bucket_size(i)]++;
    }
    for (std::map<size_t, size_t>::const_iterator it = bucket_hist.begin();
         it != bucket_hist.end(); ++it)
    {
        //fprintf(stderr, "size %zu buckets: %zu\n", (*it).first, (*it).second);
    }

}


//produce index consistent with ordering [rname, pos]
size_t SamOrder::samline_position_align(char const* samline)
{

    char contig[1024];
    size_t ones_based_pos;
    sscanf(samline, "%*[^\t]\t%*u\t%[^\t]\t%zu", contig, &ones_based_pos);

    size_t zero_based_pos = ones_based_pos == 0 ? 0 : ones_based_pos - 1;

    CONTIG_OFFSETS::const_iterator dummy;
    return flattened_position_aux(contig, zero_based_pos, this->contig_offsets, &dummy);
}


// For all records, ordering is [rname, pos]. Unmapped reads come at
// the end. For those alone, fragment_id field is sorted. This is
// achieved by adding fragment_id to the index value only when.  contig_offsets by
// convention assigns the offset of '*' (unmapped) to the last
// position in the metacontig.
size_t SamOrder::samline_projected_position_align(char const* samline)
{

    char * contig;
    char * cigar;

    size_t ones_based_pos;
    size_t flat_pos;
    CONTIG_OFFSETS::const_iterator dummy;

    sscanf(samline, "%*[^\t]\t%*u\t%a[^\t]\t%zu\t%*i\t%as", 
           &contig, &ones_based_pos, &cigar);

    size_t zero_based_pos = ones_based_pos == 0 ? 0 : ones_based_pos - 1;

    PROJECTIONS::const_iterator proj_iter = this->projections.find(contig);
    if (proj_iter == this->projections.end())
    {
        // don't do any projection
        flat_pos = flattened_position_aux(contig, zero_based_pos, 
                                          this->contig_offsets, &dummy);
    }
    else
    {
        SequenceProjection const& proj = (*proj_iter).second;
        flat_pos = flattened_position_aux(proj.target_dna.c_str(), 
                                          ExpandedStartPos(proj, zero_based_pos, cigar),
                                          this->contig_offsets, &dummy);
    }

    // now, add in the sub-ordering if contig is 
    if (strcmp(contig, "*") == 0)
    {
        flat_pos += (this->*(this->parse_fragment_id))(samline);
    }
    
    delete contig;
    delete cigar;
    return flat_pos;
}


//compute the fragment id assuming a numeric read id format
size_t SamOrder::samline_fragment_id(char const* samline)
{
    size_t fid = (this->*(this->parse_fragment_id))(samline);
    if (fid == QNAME_FORMAT_ERROR)
    {
        fprintf(stderr, "Error: samline_fragment_id: bad qname format in this line\n%s\n\n",
                samline);
        exit(1);
    }
    return fid;
}


less_seq_projection::less_seq_projection(SamOrder const* _so) :
    sam_order(_so) { }


bool less_seq_projection::operator()(SequenceProjection const& a,
                                     SequenceProjection const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;
    size_t a_pos = flattened_position_aux(a.target_dna.c_str(), a.target_start_pos(), 
                                          this->sam_order->contig_offsets, &dummy);

    size_t b_pos = flattened_position_aux(b.target_dna.c_str(), b.target_start_pos(), 
                                          this->sam_order->contig_offsets, &dummy);

    return 
        (a_pos < b_pos
         || (a_pos == b_pos
             && (a.target_end_pos() < b.target_end_pos()
                 || (a.target_end_pos() == b.target_end_pos()
                     && (strcmp(a.source_dna.c_str(), b.source_dna.c_str()) < 0)))));
}


//compute the minimum collapsed start position between the alignment
//and the guide.
size_t SamOrder::samline_position_min_align_guide(char const* samline)
{
    char guide_chrom_left[32];
    char guide_chrom_right[32];
    int read_num_left, read_num_right;
    SamFlag flag;
    size_t guide_pos_left, guide_pos_right;
    char align_chrom[32];
    size_t align_pos;

    size_t flag_raw;

    int nfields = 
        sscanf(samline, 
               "%*[^:]:" // id field.  ignored
               "read%i:%[^:]:%*c:%zu:%*[^:]:%*u:"  //read chunk
               "read%i:%[^:]:%*c:%zu:%*[^:]:%*u:" //second read chunk
               "fragment_size:%*i\t"
               "%zu\t%[^\t]\t%zu", //part of alignment
               &read_num_left, guide_chrom_left, &guide_pos_left, 
               &read_num_right, guide_chrom_right, &guide_pos_right,
               &flag_raw, align_chrom, &align_pos);
    
    if (nfields != 9)
    {
        fprintf(stderr, "SamOrder::samline_position_min_align_guide:\n"
                "Bad ID format for this sort type:\n%s\n", samline);
        exit(1);
    }

    flag.set_raw(flag_raw);

    int read_num = flag.first_fragment_in_template ? 1 : 2;

    char * guide_chrom = read_num == read_num_left ? guide_chrom_left : guide_chrom_right;
    size_t guide_pos = read_num == read_num_left ? guide_pos_left : guide_pos_right;

    CONTIG_OFFSETS::const_iterator guide_iter = this->contig_offsets.find(guide_chrom);

    if (guide_iter == this->contig_offsets.end())
    {
        fprintf(stderr, "This samline has its guide contig (%s) "
                "that is not in SAM contig index\n\n%s\n", 
                guide_chrom, samline);
        exit(1);
    }
    size_t guide_index = (*guide_iter).second + guide_pos;

    CONTIG_OFFSETS::const_iterator dummy;
    size_t align_index = 
        flattened_position_aux(align_chrom, align_pos, this->contig_offsets, &dummy);

    return std::min(guide_index, align_index);

}


// NON-member functions
//check that the genome contig order implied by the combination of
//transcript projections and transcript contig order is the same as
//that given by genome_contig_order
void GenerateProjectionHeader(CONTIG_OFFSETS const& genome_contig_order,
                              std::set<SequenceProjection> const& tx_to_genome,
                              FILE * out_sam_header_fh)
{
    //1. Compute a flattened genome coordinate for each transcript
    std::multimap<size_t, char const*> flat_genome_start;
    std::multimap<size_t, char const*>::const_iterator tit;

    std::map<char const*, size_t> lengths;

    std::set<SequenceProjection>::const_iterator pit;
    OFFSETS_ITER oit;
    for (pit = tx_to_genome.begin(); pit != tx_to_genome.end(); ++pit)
    {
        char const* tx_contig = (*pit).source_dna.c_str();
        char const* tx_proj_contig = (*pit).target_dna.c_str();
        oit = genome_contig_order.find(tx_proj_contig);
        if (oit == genome_contig_order.end())
        {
            fprintf(stderr, "Error: GenerateProjectionHeader: didn't find genome contig %s "
                    "that was found in transcript-to-genome projection from gtf file\n",
                    tx_proj_contig);
            exit(1);
        }
        size_t genome_contig_offset = (*oit).second;
        size_t tx_local_offset = (*pit).transformation.empty() ? 0 : (*pit).transformation[0].jump_length;
        size_t tx_length = (*pit).total_block_length;

        flat_genome_start.insert(std::make_pair(genome_contig_offset + tx_local_offset, tx_contig));
        lengths.insert(std::make_pair(tx_contig, tx_length));
    }

    fprintf(out_sam_header_fh, "@HD\tVN:samutil\tSO:coordinate\n");
    for (tit = flat_genome_start.begin(); tit != flat_genome_start.end(); ++tit)
    {
        fprintf(out_sam_header_fh, "@SQ\tSN:%s\tLN:%Zu\n", (*tit).second, lengths[(*tit).second]);
    }
}



//compute the fragment id assuming a numeric read id format
size_t SamOrder::parse_fragment_id_numeric(char const* qname)
{

    size_t fragment_id;
    int nfields_read = sscanf(qname, "%zu", &fragment_id);
    if (nfields_read != 1)
    {
        return QNAME_FORMAT_ERROR;
    }

    return fragment_id;
}


//compute the flattened coordinate start position.  
//read id that starts out with illumina format:
//@FLOWCELL:LANE:TILE:XPOS:YPOS
// The layout, with high bits to the left, is the same
// ffffllll tttttttt ttttxxxx xxxxxxxx xxxxxxxx xxyyyyyy yyyyyyyy yyyyyyyy
struct IlluminaID
{
    unsigned int ypos : 22;
    unsigned int xpos : 22;
    unsigned int tile : 12;
    unsigned int lane : 4;
    unsigned int flowcell : 4;
    size_t get_raw() const;
};

size_t IlluminaID::get_raw() const
{
    size_t raw = static_cast<size_t>
        (this->ypos % (1L<<22)
         | (this->xpos % (1L<<22))<<22
         | (this->tile % (1L<<12))<<44
         | (this->lane % (1L<<4))<<56
         | (this->flowcell % (1L<<4))<<60);

    return raw;
}
 // buggy test version

/*
size_t IlluminaID::get_raw() const
{
    size_t raw = static_cast<size_t>
        (this->ypos % (1L<<22)
         | (this->xpos % (1L<<22))<<22
         | (this->tile % (1L<<12))<<44
         | (this->lane % (1L<<4))<<48
         | (this->flowcell % (1L<<4))<<52);

    return raw;
}
*/


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
uint SamOrder::flowcell_hash_value(char const* flowcell)
{
    std::map<char const*, uint, less_char_ptr>::iterator fit = 
        this->flowcell_hash.find(flowcell);

    if (fit == this->flowcell_hash.end())
    {
        // need to do mutex here...
        omp_set_lock(&this->flowcell_hash_write_lock);
        size_t n = this->flowcell_hash.size();
        char * fcopy = new char[strlen(flowcell) + 1];
        strcpy(fcopy, flowcell);
        fit = this->flowcell_hash.insert(std::make_pair(fcopy, n + 1)).first;
        omp_unset_lock(&this->flowcell_hash_write_lock);
    }
    return (*fit).second;
}


// in case of format error, returns SIZE_MAX
size_t SamOrder::parse_fragment_id_illumina(char const* qname)
{

    //char flowcell[256];
    unsigned int lane, tile, xpos, ypos;

    char flowcell[256];

    IlluminaID iid;
    int end_pos;

    int nfields_read = sscanf(qname, "%[^:]:%u:%u:%u:%u%n", 
                              flowcell, &lane, &tile, &xpos, &ypos,
                              &end_pos);
    // if there is extra stuff on qname after ypos,
    // it must be one of (<spaces>, '#', '/', or '\0')
    // anything else generates an error
    if (nfields_read != 5 
        || (
            ! isspace(qname[end_pos])
            && qname[end_pos] != '#' 
            && qname[end_pos] != '/'
            && qname[end_pos] != '\0')
        )
    {
        return QNAME_FORMAT_ERROR;
    }

    iid.flowcell = this->flowcell_hash_value(flowcell);
    iid.ypos = ypos;
    iid.xpos = xpos;
    iid.tile = tile;
    iid.lane = lane;

    return iid.get_raw();
}


/*
  Parses Casava 1.8 read id format.  This routine ignores
  instrument-name, run ID, flowcell ID, and is therefore unsuitable for
  running data that comes from different flowcells etc.

  I will remedy this in the future.

  @<instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos> \
  <read number>:<is filtered>:<control number>:<barcode sequence>

  The spec is vague w.r.t header format.  But, I will assume that the statement:

  The first line is prefixed by the “@” symbol and contains the read
  name. These names are parsed until the first encountered
  whitespace. Due to this behavior, adding additional tags to the header
  line is not problematic for extant FASTQ parsers.

  means that a valid first fastq line may or may not contain a space
  followed by extra characters.

  So, another valid format could be:


  @<instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos>

*/
size_t SamOrder::parse_fragment_id_casava_1_8(char const* qname)
{

    unsigned int lane, tile, xpos, ypos;

    IlluminaID iid;

    int end_pos;
    int hash_end_pos;

    int nfields_read = sscanf(qname, "%*[^:]:%*[^:]:%*[^:]%n:%u:%u:%u:%u%n", 
                              &hash_end_pos, &lane, &tile, &xpos, &ypos, &end_pos);
    if (nfields_read != 4 ||
        (! isspace(qname[end_pos]) && qname[end_pos] != '\0'))
    {
        return QNAME_FORMAT_ERROR;
    }

    char hashable_id[256];
    strncpy(hashable_id, qname, hash_end_pos);
    hashable_id[hash_end_pos] = '\0';

    iid.flowcell = this->flowcell_hash_value(hashable_id);
    iid.lane = lane;
    iid.tile = tile;
    iid.xpos = xpos;
    iid.ypos = ypos;

    return iid.get_raw();
}



size_t SamOrder::parse_fragment_id_zero(char const* /* qname */)
{
    return 0;
}
