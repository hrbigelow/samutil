#include "sam_order.h"
#include "sam_helper.h"

#include <set>
#include <zlib.h>

//!!! This is muddied in that it is partially up to the user, partly a nuissance.
SamOrder::SamOrder(SAM_ORDER _so, char const* _sort_type) : 
    order(_so)
{
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
    else if (strcmp(this->sort_type, "READ_ID_FLAG") == 0)
    {
        this->sam_index = &SamOrder::samline_read_id;
    }
    else if (strcmp(this->sort_type, "FRAGMENT_ID") == 0)
    {
        this->sam_index = &SamOrder::samline_fragment_id;
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


SamOrder::SamOrder() : order(SAM_RID)
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

    // for (size_t k = 0; k != contig_offsets.size(); ++k)
    // {
    //     delete keys[k];
    // }
    // delete keys;
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


//Call to initialize the proper READ ID parsing format
SAM_QNAME_FORMAT SamOrder::InitFromFastqFile(char const* file)
{
    
    char id[1024];

    gzFile fh = gzopen(file, "r");
    gzread(fh, id, 1023);
    gzclose(fh);

    // a fastq file has a single first line with '@' followed by the ID.
    // that's why there is a '1' here and in the next statement.
    char * nl = strchr(id, '\n');
    nl == NULL ? id[1] = '\0' : *nl = '\0';

    return InitFromID(id + 1);

}


SAM_QNAME_FORMAT SamOrder::InitFromID(char const* id)
{

    int dum[5];
    char dum_string[30];

    SAM_QNAME_FORMAT qname_format;

    if (strlen(id) == 0)
    {
        this->parse_fragment_id = &parse_fragment_id_zero;
        qname_format = SAM_NON_INTERPRETED;
    }
    else if (sscanf(id, "%i", dum) == 1)
    {
        this->parse_fragment_id = &parse_fragment_id_numeric;
        qname_format = SAM_NUMERIC;
    }
    else if (sscanf(id, "%[^:]:%i:%i:%i:%i", dum_string, dum, dum + 1, dum + 2, dum + 3) == 5)
    {
        this->parse_fragment_id = &parse_fragment_id_illumina;
        qname_format = SAM_ILLUMINA;
    }

    else if (sscanf(id, "%*[^:]:%*[^:]:%[^:]:%i:%i:%i:%i", 
                    dum_string, dum, dum + 1, dum + 2, dum + 3) == 5)

    {
        this->parse_fragment_id = &parse_fragment_id_casava_1_8;
        qname_format = SAM_CASAVA18;
    }
    else
    {
        fprintf(stderr, "Error: SamOrder: Don't have a parser for this qname format: %s\n",
                id);
        exit(1);
    }

    return qname_format;
}


void SamOrder::InitFromChoice(SAM_QNAME_FORMAT qname_format)
{
    switch(qname_format)
    {
    case SAM_NON_INTERPRETED: this->parse_fragment_id = &parse_fragment_id_zero; break;
    case SAM_NUMERIC: this->parse_fragment_id = &parse_fragment_id_numeric; break;
    case SAM_ILLUMINA: this->parse_fragment_id = &parse_fragment_id_illumina; break;
    case SAM_CASAVA18: this->parse_fragment_id = &parse_fragment_id_casava_1_8; break;
    default:
        fprintf(stderr, "Error: InitFromChoice: Don't know this qname format\n");
        exit(1);
    }
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
        fprintf(stderr, "flattened_position: error: rname %s (at %Zu) "
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
           (strcmp(a.cigar, b.cigar) < 0 ||
            (strcmp(a.cigar, b.cigar) == 0 &&
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
        strcmp(a.cigar, b.cigar) == 0 &&
        flattened_position_mate(&a, &dummy) == flattened_position_mate(&b, &dummy) &&
        a.tags.alignment_space == b.tags.alignment_space;

        // flattened_position(&a, &dummy) == flattened_position(&b, &dummy) &&
        // a.this_fragment_on_pos_strand() == b.this_fragment_on_pos_strand() &&
        // strcmp(a.cigar, b.cigar) == 0 &&
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
           (strcmp(a.cigar, b.cigar) < 0))));

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
void SamOrder::AddHeaderContigStats(FILE * sam_fh)
{
    size_t contig_offset = 0;

    char line[4096 + 1];

    while (! feof(sam_fh))
    {
        fgets(line, 4096, sam_fh);

        if (line[0] != '@')
        {
            rewind(sam_fh);
            break;
        }
        //if there is a SQ tag, parse and record
        else if (strncmp(line, "@SQ", 3) == 0)
        {
            //this is a SQ line
            char contig_name[1000];
            
            size_t contig_length;
            sscanf(line, "@SQ\tSN:%s\tLN:%zu", contig_name, &contig_length);
            this->contig_lengths[std::string(contig_name)] = contig_length;

            char * contig_name_copy = new char[strlen(contig_name) + 1];
            strcpy(contig_name_copy, contig_name);

            this->contig_offsets[contig_name_copy] = contig_offset;
            contig_offset += contig_length;
        }
    }

    this->contig_lengths[std::string("*")] = 0;

    
    char * contig_name_copy = new char[2];
    strcpy(contig_name_copy, "*");

    this->contig_offsets[contig_name_copy] = 0;

    if (contig_offset == 0)
    {
        fprintf(stderr, "No SQ lines in header, or no SN / LN fields defining contigs\n");
        exit(1);
    }

    this->contig_offsets.rehash(this->contig_offsets.size() * 10);

    fprintf(stderr, "size()=%zu, bucket_count()=%zu, load_factor()=%f\n", 
            this->contig_offsets.size(),
            this->contig_offsets.bucket_count(),
            this->contig_offsets.load_factor());

    std::map<size_t, size_t> bucket_hist;
    for (size_t i = 0; i != this->contig_offsets.bucket_count(); ++i)
    {
        bucket_hist[this->contig_offsets.bucket_size(i)]++;
    }
    for (std::map<size_t, size_t>::const_iterator it = bucket_hist.begin();
         it != bucket_hist.end(); ++it)
    {
        fprintf(stderr, "size %zu buckets: %zu\n",
                (*it).first, (*it).second);
    }

}


//produce index consistent with ordering [rname, pos]
size_t SamOrder::samline_position_align(char const* samline) const
{

    char contig[1024];
    size_t position;
    sscanf(samline, "%*s\t%*u\t%s\t%zu", contig, &position);

    CONTIG_OFFSETS::const_iterator dummy;
    return flattened_position_aux(contig, position, this->contig_offsets, &dummy);
}


//compute the fragment id assuming a numeric read id format
size_t SamOrder::samline_read_id(char const* samline) const
{

    char qname[1024];
    SamFlag flag;

    // qname flag rname pos mapq cigar rnext pnext tlen seq qual tags...
    int nfields_read = sscanf(samline, "%s\t%zu", qname, &flag.raw);
    if (nfields_read != 2)
    {
        fprintf(stderr, "samline_read_id_flag: bad format for read-id sorting.\n"
                "Encountered data line:\n%s\n",
                samline);
        exit(1);
    }

    size_t read_id = this->parse_fragment_id(qname);

    read_id = read_id<<1;
    int read_num = flag.first_fragment_in_template ? 0 : 1;
    read_id += read_num;

    return read_id;
}



//compute the fragment id assuming a numeric read id format
size_t SamOrder::samline_fragment_id(char const* samline) const
{
    return this->parse_fragment_id(samline);
}


//compute the minimum collapsed start position between the alignment
//and the guide.
size_t SamOrder::samline_position_min_align_guide(char const* samline) const
{
    char guide_chrom_left[32];
    char guide_chrom_right[32];
    int read_num_left, read_num_right;
    SamFlag flag;
    size_t guide_pos_left, guide_pos_right;
    char align_chrom[32];
    size_t align_pos;

    sscanf(samline, 
           "%*[^:]:" // id field.  ignored
           "read%i:%[^:]:%*c:%zu:%*[^:]:%*u"  //read chunk
           "read%i:%[^:]:%*c:%zu:%*[^:]:%*u\t" //second read chunk
           "%zu\t%s\t%zu", //part of alignment
           &read_num_left, guide_chrom_left, &guide_pos_left, 
           &read_num_right, guide_chrom_right, &guide_pos_right,
           &flag.raw, align_chrom, &align_pos);
    
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
size_t parse_fragment_id_numeric(char const* qname)
{

    size_t fragment_id;
    int nfields_read = sscanf(qname, "%zu", &fragment_id);
    if (nfields_read != 1)
    {
        fprintf(stderr, "parse_fragment_id_numeric: "
                "Error: qname format is not numeric: %s\n", qname);
        exit(1);
    }

    return fragment_id;
}



//compute the flattened coordinate start position.  
//read id that starts out with illumina format:
//@FLOWCELL:LANE:TILE:XPOS:YPOS
union IlluminaID
{
    size_t data;
    struct
    {
        unsigned int ypos : 16;
        unsigned int xpos : 16;
        unsigned int tile : 8;
        unsigned int lane : 8;
        unsigned int flowcell : 16;
    } f;
};


//update record of flow cells.
//contains a static member!
//is this considered a memory leak?
uint flowcell_hash_value(char const* flowcell)
{
    static std::map<char const*, uint, less_char_ptr> flowcell_hash;

    std::map<char const*, uint, less_char_ptr>::iterator fit = flowcell_hash.find(flowcell);
    if (fit == flowcell_hash.end())
    {
        size_t n = flowcell_hash.size();
        char * fcopy = new char[strlen(flowcell) + 1];
        strcpy(fcopy, flowcell);
        fit = flowcell_hash.insert(std::make_pair(fcopy, n + 1)).first;
    }
    return (*fit).second;
}


size_t parse_fragment_id_illumina(char const* qname)
{

    //char flowcell[256];
    uint lane, tile, xpos, ypos;

    char flowcell[256];

    IlluminaID iid;

    int nfields_read = sscanf(qname, "%[^:]:%i:%i:%i:%i", flowcell, &lane, &tile, &xpos, &ypos);
    if (nfields_read != 5)
    {
        fprintf(stderr, "parse_fragment_id_illumina: bad format for read-id sorting\n");
        exit(1);
    }

    iid.f.flowcell = flowcell_hash_value(flowcell);
    iid.f.lane = lane;
    iid.f.tile = tile;
    iid.f.xpos = xpos;
    iid.f.ypos = ypos;

    return iid.data;
}


/*
Parses Casava 1.8 read id format.  This routine ignores
instrument-name, run ID, flowcell ID, and is therefore unsuitable for
running data that comes from different flowcells etc.

I will remedy this in the future.

@ <instrument-name>:<run ID>:<flowcell ID>:<lane-number>:<tile-number>:<x-pos>:<y-pos> \
<read number>:<is filtered>:<control number>:<barcode sequence>
*/
size_t parse_fragment_id_casava_1_8(char const* qname)
{

    uint lane, tile, xpos, ypos;

    char flowcell[256];

    IlluminaID iid;

    int nfields_read = sscanf(qname, "%*[^:]:%*[^:]:%[^:]:%i:%i:%i:%i", flowcell, &lane, &tile, &xpos, &ypos);
    if (nfields_read != 5)
    {
        fprintf(stderr, "parse_fragment_id_casava_1_8: bad format for read-id sorting\n");
        exit(1);
    }

    iid.f.flowcell = flowcell_hash_value(flowcell);
    iid.f.lane = lane;
    iid.f.tile = tile;
    iid.f.xpos = xpos;
    iid.f.ypos = ypos;

    return iid.data;
}



size_t parse_fragment_id_zero(char const* /* qname */)
{
    return 0;
}
