#include "sam_order.h"
#include "sam_helper.h"

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
        this->less = &SamOrder::less_fid_position; 
        this->equal = &SamOrder::equal_fid_position; 
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
        this->sam_index = &SamOrder::samline_read_id_flag;
    }
    else
    {
        this->sam_index = NULL;
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
    char const** keys = new char const*[contig_offsets.size()];
    size_t k;
    CONTIG_OFFSETS::iterator ci;
    for (ci = contig_offsets.begin(), k = 0; ci != contig_offsets.end(); ++ci, ++k)
    {
        keys[k] = (*ci).first;
    }
    contig_offsets.clear();

    for (size_t k = 0; k != contig_offsets.size(); ++k)
    {
        delete keys[k];
    }
    delete keys;
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
inline int fid_aux(SamLine const& a, SamLine const& b)
{
    int qname_cmp;
    if (! (qname_cmp = static_cast<int>(a.qid) - static_cast<int>(b.qid)))
    {
        qname_cmp = strcmp(a.qname, b.qname);
    }
    return qname_cmp;
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
    
    return flattened_position_aux(a->mate_ref_name(), a->mpos,
                                  this->contig_offsets,
                                  contig_iter);
}




//for distinguishing SamLines by alignment position,
//[rname, pos, query_strand, cigar, mrnm, mpos]
bool SamOrder::less_position(SamLine const& a, SamLine const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;
    int64_t flattened_cmp = 
        static_cast<int64_t>(this->flattened_position(&a, &dummy))
        - static_cast<int64_t>(this->flattened_position(&b, &dummy));

    return flattened_cmp < 0 || 
        (flattened_cmp == 0 && 
         (a.query_on_pos_strand() < b.query_on_pos_strand() ||
          (a.query_on_pos_strand() == b.query_on_pos_strand() &&
           (strcmp(a.cigar, b.cigar) < 0 ||
            (strcmp(a.cigar, b.cigar) == 0 &&
             (static_cast<int64_t>(this->flattened_position_mate(&a, &dummy))
              - static_cast<int64_t>(this->flattened_position_mate(&b, &dummy)) < 0))))));
    
}

//[rname, pos, query_strand, cigar, mrnm, mpos]
bool SamOrder::equal_position(SamLine const& a, SamLine const& b) const
{

    CONTIG_OFFSETS::const_iterator dummy;
    int64_t flattened_cmp = 
        static_cast<int64_t>(flattened_position(&a, &dummy)) 
        - static_cast<int64_t>(flattened_position(&b, &dummy));

    return 
        flattened_cmp == 0 &&
        a.query_on_pos_strand() == b.query_on_pos_strand() &&
        strcmp(a.cigar, b.cigar) == 0 &&
        flattened_position_mate(&a, &dummy) == flattened_position_mate(&b, &dummy);
}


//for distinguishing SamLines by fragment alignment position
//[rname, pos, query_strand, cigar, mrnm, mpos]
bool SamOrder::less_fposition(SamLine const& a, SamLine const& b) const
{
    CONTIG_OFFSETS::const_iterator dummy;

    size_t a_fragment_pos = a.mapped_in_proper_pair()
        ? std::min(this->flattened_position(&a, &dummy),
                   this->flattened_position_mate(&a, &dummy))
        : this->flattened_position(&a, &dummy);

    size_t b_fragment_pos = b.mapped_in_proper_pair()
        ? std::min(this->flattened_position(&b, &dummy),
                   this->flattened_position_mate(&b, &dummy))
        : this->flattened_position(&b, &dummy);

    int64_t flattened_cmp = static_cast<int64_t>(a_fragment_pos)
        - static_cast<int64_t>(b_fragment_pos);

    return flattened_cmp < 0 || 
        (flattened_cmp == 0 && 
         (a.query_on_pos_strand() < b.query_on_pos_strand() ||
          (a.query_on_pos_strand() == b.query_on_pos_strand() &&
           (strcmp(a.cigar, b.cigar) < 0))));

}

//consider the 'read id' to be composed of its qname and pair identity flag
//[qname, pair]
bool SamOrder::less_rid(SamLine const& a, SamLine const& b) const
{
    int qname_cmp = fid_aux(a, b);
    return qname_cmp < 0 
        || (qname_cmp == 0
            && (a.first_read_in_pair() < b.first_read_in_pair()));
}


//[qname, pair]
bool SamOrder::equal_rid(SamLine const& a, SamLine const& b) const
{
    int qname_cmp = fid_aux(a, b);
    return qname_cmp == 0
        && (a.first_read_in_pair() == b.first_read_in_pair());
}


//[qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
bool SamOrder::less_rid_position(SamLine const& a, SamLine const& b) const
{
    return less_rid(a, b) ||
        (equal_rid(a, b) && less_position(a, b));
}

//[qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
bool SamOrder::equal_rid_position(SamLine const& a, SamLine const& b) const
{
    return equal_rid(a, b) && equal_position(a, b);
}


//[rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
bool SamOrder::less_position_rid(SamLine const& a, SamLine const& b) const
{
    return less_position(a, b) ||
        (equal_position(a, b) && less_rid(a, b));
}

//[rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
bool SamOrder::equal_position_rid(SamLine const& a, SamLine const& b) const
{
    return equal_position(a, b) && equal_rid(a, b);
}




bool SamOrder::less_fid(SamLine const& a, SamLine const& b) const
{
    return fid_aux(a, b) < 0;
}


bool SamOrder::equal_fid(SamLine const& a, SamLine const& b) const
{
    return fid_aux(a, b) == 0;
}

bool SamOrder::less_fid_position(SamLine const& a, SamLine const& b) const
{
    return less_fid(a, b)
        || (equal_fid(a, b)
            && less_position(a, b));
}

bool SamOrder::equal_fid_position(SamLine const& a, SamLine const& b) const
{
    return equal_fid(a, b) && equal_position(a, b);
}




//compute contig lengths from sam_fh.
//rewinds sam_fh to beginning
void SamOrder::AddHeaderContigStats(FILE * sam_fh)
{
    SamLine * samline;
    size_t contig_offset = 0;
    bool allow_absent_seq_qual = true;
    while ((samline = new SamLine(sam_fh, allow_absent_seq_qual)))
    {
        if (samline->parse_flag != HEADER)
        {
            //SamLine::line does not have a newline at the end.  But we want to count it.
            delete samline;
            rewind(sam_fh);
            break;
            // long last_line_length = strlen(samline->line) + 1;
            // fseek(sam_fh, - last_line_length, std::ios::cur);
            // break;
        }
        //if there is a SQ tag, parse and record
        else if (strncmp(samline->tag_string, "SQ", 2) == 0)
        {
            //this is a SQ line
            char contig_name[1000];
            
            size_t contig_length;
            sscanf(samline->tag_string, "SQ\tSN:%s\tLN:%zu", contig_name, &contig_length);
            this->contig_lengths[std::string(contig_name)] = contig_length;

            char * contig_name_copy = new char[strlen(contig_name) + 1];
            strcpy(contig_name_copy, contig_name);

            this->contig_offsets[contig_name_copy] = contig_offset;
            contig_offset += contig_length;
        }
        delete samline;
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

}


//compute the flattened coordinate start position

//produce index consistent with ordering [rname, pos]
size_t SamOrder::samline_position_align(char const* samline) const
{

    char contig[1024];
    size_t position;
    sscanf(samline, "%*s\t%*u\t%s\t%zu", contig, &position);

    CONTIG_OFFSETS::const_iterator dummy;
    return flattened_position_aux(contig, position, this->contig_offsets, &dummy);
}


//compute the flattened coordinate start position.  only works for a read id that
//starts out numerically
size_t SamOrder::samline_read_id_flag(char const* samline) const
{

    char qname[1024];
    size_t flag;
    size_t qnumber;

    // qname flag rname pos mapq cigar mrnm mpos isize seq qual tags...
    int nfields_read = sscanf(samline, "%s\t%zu", qname, &flag);
    if (nfields_read != 2)
    {
        fprintf(stderr, "samline_read_id_flag: bad format for read-id sorting\n");
        exit(1);
    }
    nfields_read = sscanf(qname, "%zu", &qnumber);
    if (nfields_read != 1)
    {
        fprintf(stderr, "samline_read_id_flag: bad format for read-id sorting\n");
        exit(1);
    }

    qnumber = qnumber<<1;
    int read_num = ((flag & SamFlags::FIRST_READ_IN_PAIR) != 0) ? 0 : 1;
    qnumber += read_num;

    return qnumber;
}


//index consistent with sort order: [qname, pair, rname, pos].
//This indexing stains the capacity of size_t.  We want to allow
//2^32 bits for qname, pair, and another 2^32 for rname, pos.
//It is not likely but possible that in the future, either of these limits will be exceeded.

/*
size_t index_sam_fid_position(char const* samline, CONTIG_OFFSETS const& contig_offsets)
{
    // qname flag rname pos mapq cigar mrnm mpos isize seq qual tags...
    char qname[1024];
    size_t flag;
    size_t qid;

    int nfields_read = sscanf(samline, "%zu%*s\t%zu\t%s\t%zu", &qid, &flag, rname, &pos);
    if (nfields_read != 4)
    {
        fprintf(stderr, "samline_read_id_flag: bad format for rid_position sorting\n");
        exit(1);
    }


    qid = qid<<1;

    if (qid >= 2<<32)
    {
        fprintf(stderr, "Error: index_sam_rid_position: qid > 2^32.  cannot index more than 2^32 fragments.\n");
        exit(1);
    }

    int read_num = ((flag & SamFlags::FIRST_READ_IN_PAIR) != 0) ? 0 : 1;
    qid += read_num;
    qid = qid<<32;

    CONTIG_OFFSETS::const_iterator contig_iter = contig_offsets.find(rname);

    if (contig_iter == contig_offsets.end())
    {
        fprintf(stderr, "This samline has its guide contig (%s) "
                "that is not in SAM contig index\n\n%s\n", 
                rname, samline);
        exit(1);
    }
    size_t contig_index = (*contig_iter).second + pos;

    if (qid >= 2<<32 || contig_index >= 2<<32)
    {
        fprintf(stderr, "Error: index_sam_rid_position: in indexing rid %zu, rname %s, pos %zu, (position index %zu)"
                "capacity exceeded.  Cannot index more than 2^32 fragments or alignment positions.\n",
                qid, rname, pos, contig_index);
        exit(1);
    }
    return qid & contig_index;
}
*/

//compute the minimum collapsed start position between the alignment
//and the guide.
size_t SamOrder::samline_position_min_align_guide(char const* samline) const
{
    char guide_chrom_left[32];
    char guide_chrom_right[32];
    int read_num_left, read_num_right;
    size_t flag;
    size_t guide_pos_left, guide_pos_right;
    char align_chrom[32];
    size_t align_pos;

    sscanf(samline, 
           "%*[^:]:" // id field.  ignored
           "read%i:%[^:]:%*c:%zu:%*[^:]:"  //read chunk
           "read%i:%[^:]:%*c:%zu:%*s\t" //second read chunk
           "%zu\t%s\t%zu", //part of alignment
           &read_num_left, guide_chrom_left, &guide_pos_left, 
           &read_num_right, guide_chrom_right, &guide_pos_right,
           &flag, align_chrom, &align_pos);
    
    int read_num = ((flag & SamFlags::FIRST_READ_IN_PAIR) != 0) ? 1 : 2;

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
