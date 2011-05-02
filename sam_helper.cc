#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "sam_helper.h"
#include "cigar_ops.h"
#include "file_utils.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <fstream>

namespace SamFlags
{
    int const PAIRED_IN_SEQUENCING = 0x0001;
    int const MAPPED_IN_PROPER_PAIR = 0x0002;
    int const QUERY_UNMAPPED = 0x0004;
    int const MATE_UNMAPPED = 0x0008;
    int const QUERY_ON_NEG_STRAND = 0x0010; //set when negative stranded
    int const MATE_ON_NEG_STRAND = 0x0020; //set when mate is negative stranded
    int const FIRST_READ_IN_PAIR = 0x0040;
    int const SECOND_READ_IN_PAIR = 0x0080;
    int const ALIGNMENT_NOT_PRIMARY = 0x0100;
    int const FAILED_QUALITY_CHECK = 0x0200;
    int const PCR_OR_OPTICAL_DUPLICATE = 0x0400;
};


//serialize the data in a line.  string data will be packed in the
//'line' variable, with null's terminating each one.  non-string data
//will be stored separately

bool SamLine::numeric_start_fragment_ids = false;

SamLine::SamLine(SAM_PARSE _parse_flag,
                 char const* _qname, int _flag, 
                 char const* _rname, size_t _pos,
                 size_t _mapq, char const* _cigar,
                 char const* _mrnm, size_t _mpos,
                 int _isize, char const* _seq,
                 char const* _qual,
                 char const* _tag_string) :
    parse_flag(_parse_flag), flag(_flag), pos(_pos), 
    mapq(_mapq), mpos(_mpos), 
    isize(_isize), extra(NULL), extra_tag(NULL)
{
    
    assert(strlen(_qual) < 10000);
    assert(strlen(_seq) < 10000);

    char const* src[] = { 
        _qname, _rname, _cigar, _mrnm, _seq, _qual, _tag_string
    };
    size_t len[7];

    this->bytes_in_line = 0;
    for (size_t f = 0; f != 7; ++f)
    {
        len[f] = strlen(src[f]);
        this->bytes_in_line += len[f];
    }
    this->bytes_in_line += 7;

    char ** trg[] = { 
        &this->qname, &this->rname, &this->cigar, 
        &this->mrnm, &this->seq, &this->qual,
        &this->tag_string
    };

    this->line = new char[this->bytes_in_line + 1];

    this->line[0] = '\0';
    char * dest = this->line;
    for (size_t f = 0; f != 7; ++f)
    {
        strcpy(dest, src[f]);
        *trg[f] = dest;
        dest += len[f] + 1;
    }
    
    //now, make empty seq, qual, and tag_string fields NULL if they are empty
    if (strlen(_seq) == 0)
    {
        this->seq = NULL;
    }
    if (strlen(_qual) == 0)
    {
        this->qual = NULL;
    }
    if (strlen(_tag_string) == 0)
    {
        this->tag_string = NULL;
    }

    if (SamLine::numeric_start_fragment_ids)
    {
        int nfields_read = sscanf(this->qname, "%zu", &this->qid);
        if (nfields_read != 1)
        {
            fprintf(stderr, "SamLine: Error: numeric_start_fragment_ids set to true, but qname does not "
                    "start with integer:\n%s", this->qname);
            exit(1);
        }
    }
}


//for distinguishing SamLines by alignment position,
//[rname, pos, query_strand, cigar, mrnm, mpos]
bool SamLine::less_position(SamLine const& b) const
{
    int refcmp = strcmp(this->rname, b.rname);

    return refcmp < 0 || 
        (refcmp == 0 && 
         (this->pos < b.pos || 
          (this->pos == b.pos &&
           (this->query_on_pos_strand() < b.query_on_pos_strand() ||
            (this->query_on_pos_strand() == b.query_on_pos_strand() &&
             (strcmp(this->cigar, b.cigar) < 0 ||
              (strcmp(this->cigar, b.cigar) == 0 &&
               (strcmp(this->mrnm, b.mrnm) < 0 ||
                (strcmp(this->mrnm, b.mrnm) == 0 &&
                 (this->mpos < b.mpos))))))))));

}

//[rname, pos, query_strand, cigar, mrnm, mpos]
bool SamLine::equal_position(SamLine const& b) const
{
    int refcmp = strcmp(this->rname, b.rname);

    return 
        refcmp == 0 &&
        this->pos == b.pos &&
        this->query_on_pos_strand() == b.query_on_pos_strand() &&
        strcmp(this->cigar, b.cigar) == 0 &&
        strcmp(this->mrnm, b.mrnm) == 0 &&
        this->mpos == b.mpos;
}

//consider the 'read id' to be composed of its qname and pair identity flag
//[qname, pair]
bool SamLine::less_rid(SamLine const& b) const
{
    int qname_cmp = 
        SamLine::numeric_start_fragment_ids 
        ? (static_cast<int>(this->qid) - static_cast<int>(b.qid)) 
        : strcmp(this->qname, b.qname);

    return qname_cmp < 0 
        || (qname_cmp == 0
            && (this->first_read_in_pair() < b.first_read_in_pair()));
}


//[qname, pair]
bool SamLine::equal_rid(SamLine const& b) const
{
    int qname_cmp = 
        SamLine::numeric_start_fragment_ids 
        ? (static_cast<int>(this->qid) - static_cast<int>(b.qid)) 
        : strcmp(this->qname, b.qname);

    return qname_cmp == 0
        && (this->first_read_in_pair() == b.first_read_in_pair());
}


//[qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
bool SamLine::less_rid_position(SamLine const& b) const
{
    return this->less_rid(b) ||
        (this->equal_rid(b) && this->less_position(b));
}

//[qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
bool SamLine::equal_rid_position(SamLine const& b) const
{
    return this->equal_rid(b) && this->equal_position(b);
}



//[rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
bool SamLine::less_position_rid(SamLine const& b) const
{
    return this->less_position(b) ||
        (this->equal_position(b) && this->less_rid(b));
}

//[rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
bool SamLine::equal_position_rid(SamLine const& b) const
{
    return this->equal_position(b) && this->equal_rid(b);
}


bool SamLine::less_fid(SamLine const& b) const
{
    int qname_cmp = 
        SamLine::numeric_start_fragment_ids 
        ? (static_cast<int>(this->qid) - static_cast<int>(b.qid)) 
        : strcmp(this->qname, b.qname);

    return qname_cmp < 0;
}


bool SamLine::equal_fid(SamLine const& b) const
{
    int qname_cmp = 
        SamLine::numeric_start_fragment_ids 
        ? (static_cast<int>(this->qid) - static_cast<int>(b.qid)) 
        : strcmp(this->qname, b.qname);

    return qname_cmp == 0;
}


bool SamLine::less_fid_position(SamLine const& b) const
{
    return this->less_fid(b)
        || (this->equal_fid(b)
            && this->less_position(b));
}


bool SamLine::equal_fid_position(SamLine const& b) const
{
    return this->equal_fid(b) && this->equal_position(b);
}



//sort by alignment location of the read itself: (refname, position, strand, cigar)

//warning: if you need to keep track of unique read pairs in which two
//read pairs will perhaps share a given identical read,
/*
bool SamLine::operator<(SamLine const& b) const
{
    int refcmp = strcmp(this->rname, b.rname);

    return refcmp < 0 || 
        (refcmp == 0 && 
         (this->pos < b.pos || 
          (this->pos == b.pos &&
           (this->query_on_pos_strand() < b.query_on_pos_strand() ||
            (this->query_on_pos_strand() == b.query_on_pos_strand() ||
             (strcmp(this->cigar, b.cigar) < 0))))));
            
}
*/



SamLine::SamLine(FILE * seqfile, bool sam_is_ones_based, bool allow_absent_seq_qual) :
    bytes_in_line(0), line(NULL), qname(NULL), qid(0), rname(NULL), 
    cigar(NULL), mrnm(NULL), seq(NULL), qual(NULL), tag_string(NULL), extra(NULL),
    extra_tag(NULL)
{
    char buf[4096 + 1];
    if (fgets(buf, 4096, seqfile) == NULL)
    {
        this->parse_flag = END_OF_FILE;
        return;
    }
    this->Init(buf, sam_is_ones_based, allow_absent_seq_qual);
}



SamLine::SamLine(char const* samline_string, 
                 bool sam_is_ones_based,
                 bool allow_absent_seq_qual) :
    qname(NULL), qid(0), rname(NULL), cigar(NULL), mrnm(NULL), 
    seq(NULL), qual(NULL), tag_string(NULL), extra(NULL), extra_tag(NULL)
{
    this->Init(samline_string, sam_is_ones_based, allow_absent_seq_qual);
}


//samline_string must be terminated by newline
void SamLine::Init(char const* samline_string,
                   bool /* sam_is_ones_based */,
                   bool allow_absent_seq_qual)
{

    this->bytes_in_line = strlen(samline_string);
    this->line = new char[this->bytes_in_line + 1];
    
    //this->line will be double-NULL terminated.  one NULL to replace
    //the newline, and one following that
    strcpy(this->line, samline_string);
    this->line[this->bytes_in_line - 1] = '\0'; //set newline to null as well
    
    if (samline_string[0] == '@')
    {
        //this is a header
        this->parse_flag = HEADER;

        this->tag_string = this->line + 1;
        return;
    }

    else
    {
        bool was_error = false;

        //there will either be an immediate newline, or a tab, then
        //some non-tab characters, then a newline
        int num_numeric_fields =
            sscanf(this->line, "%*s\t%i\t%*s\t%zu\t%zu\t%*s\t%*s\t%zu\t%i",
                   &this->flag, &this->pos, &this->mapq, &this->mpos, &this->isize);

        if (num_numeric_fields != 5)
        {
            was_error = true;
        }

        //now, replace tabs with nulls, and set addresses
        char * ptr = this->line;

        char * dummy = NULL;
        char * end_of_targets = NULL;
        char ** trg[] = { 
            &this->qname, &dummy, &this->rname, &dummy, &dummy, &this->cigar, &this->mrnm,
            &dummy, &dummy, &this->seq, &this->qual, &this->tag_string
        };

        char *** trg_ptr = trg;
        **trg_ptr = ptr;

        //this->line is NULL-terminated
        //replace tabs with NULLs in this->line
        //set trg[] members to point to each field
        while (*ptr != '\0' && *trg_ptr != &this->tag_string)
        {
            ptr = strchrnul(ptr, '\t'); //advance to end of this field
            *ptr++ = '\0'; //replace with NULL, advance to start of next field
            **++trg_ptr = ptr; //set target to point to this field
        }
        
        bool seq_absent = *trg_ptr == &this->seq;

        if (*trg_ptr != &this->tag_string &&
            ! (allow_absent_seq_qual && seq_absent))
        {
            was_error = true;
        }

        if (seq_absent)
        {
            this->seq = NULL;
            this->qual = NULL;
            this->tag_string = NULL;
        }
        else
        {
            //seq/qual are present, is tag string there?
            if (strlen(this->tag_string) == 0)
            {
                this->tag_string = NULL;
            }
        }
        

        if (! this->query_unmapped())
        {
            --this->pos;
        }
        if (! this->mate_unmapped())
        {
            --this->mpos;
        }            

        if (SamLine::numeric_start_fragment_ids)
        {
            int nfields_read = sscanf(this->qname, "%zu", &this->qid);
            if (nfields_read != 1)
            {
            fprintf(stderr, "SamLine: Error: numeric_start_fragment_ids set to true, but qname does not "
                    "start with integer:\n%s", this->qname);
            exit(1);
            }
        }

        if (was_error)
        {
            fprintf(stderr, "SamLine::load_line: encountered badly formatted line.\n");
            fprintf(stderr, "Line:\n\n%s\n", this->line);
            this->parse_flag = PARSE_ERROR;
            return;
        }
        
        this->parse_flag = DATA_LINE;
    }
}


SamLine::~SamLine()
{
    if (this->line != NULL)
    {
        delete this->line;
        this->line = NULL;
    }
    if (this->extra != NULL)
    {
        delete this->extra;
        this->extra = NULL;
    }
    if (this->extra_tag != NULL)
    {
        delete this->extra_tag;
        this->extra_tag = NULL;
    }
}



void SamLine::print(FILE * seqfile, bool flip_query_strand_flag) const
{

    if (this->parse_flag == HEADER)
    {
        fprintf(seqfile, "%s\n", this->line);
        return;
    }

    int offset = 1; // output is always ones-based in SAM file
    int flag_to_print = flip_query_strand_flag ? 
        (this->flag ^ SamFlags::QUERY_ON_NEG_STRAND) : this->flag;
    
    fprintf(seqfile, "%s\t%i\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i",
            this->qname, flag_to_print, this->rname, this->pos + offset,
            this->mapq, this->cigar, this->mrnm, this->mpos + offset,
            this->isize);

    if (this->seq != NULL)
    {
        assert(this->qual != NULL);
        fprintf(seqfile, "\t%s\t%s", this->seq, this->qual);
    }

    if (this->tag_string != NULL)
    {
        fprintf(seqfile, "\t%s", this->tag_string);
    }
    fprintf (seqfile, "\n");
}



size_t SamLine::aligned_read_length()
{
    Cigar::CIGAR_VEC cigarvec = Cigar::FromString(this->cigar, 0);
    bool use_top_coord = true;
    return Cigar::Length(cigarvec, ! use_top_coord);
}

int SamLine::alignment_score(char const* tag, int default_if_missing, bool * has_score) const
{
    int align_score;
    char * score_string = this->tag_string == NULL ? NULL : strstr(this->tag_string, tag);

    char tag_plus_format[32];
    strcpy(tag_plus_format, tag);
    strcat(tag_plus_format, ":i:%i");

    if (score_string != NULL 
        && sscanf(score_string, tag_plus_format, &align_score) == 1)
    {
        *has_score = true;
        return align_score;
    }
    else
    {
        *has_score = false;
        return default_if_missing;
    }
}


void SamLine::add_tag(char const* tag)
{
    if (this->extra_tag == NULL)
    {
        //tag_string was previously pointing into this->line
        this->extra_tag = new char[strlen(this->tag_string) + strlen(tag) + 1];
        strcpy(this->extra_tag, this->tag_string);
        this->tag_string = this->extra_tag;
    }
    else
    {
        //extra_tag already in use.  grow it.
        char * old_extra_tag = this->extra_tag;
        this->extra_tag = new char[strlen(this->tag_string) + strlen(tag) + 1];
        strcpy(this->extra_tag, old_extra_tag);
        delete old_extra_tag;
    }
    if (strlen(this->extra_tag) != 0)
    {
        strcat(this->extra_tag, "\t");
    }
    strcat(this->extra_tag, tag);
}

//true if samlines represent sequence mate pairs, regardless of
//whether they are mapped as such
bool AreSequencedMatePairs(SamLine const& a, SamLine const& b)
{
    return 
        strcmp(a.qname, b.qname) == 0
        && (a.first_read_in_pair() != b.first_read_in_pair())
        && a.paired_in_sequencing()
        && b.paired_in_sequencing();
}

//true if samlines are sequenced mate pairs and the alignments of a
//and b are coordinated as a single alignment hypothesis. This is
//ill-defined, because it is possible to construct two distinct
//fragment alignments (of the same fragment) whose individual reads
//are ambiguous mates.
bool AreMappedMatePairs(SamLine const& a, SamLine const& b)
{

    char const* amate_ref = strcmp(a.mrnm, "=") == 0 ? a.rname : a.mrnm;
    char const* bmate_ref = strcmp(b.mrnm, "=") == 0 ? b.rname : b.mrnm;

    return
        AreSequencedMatePairs(a, b)
        && a.mpos == b.pos
        && b.mpos == a.pos
        && (a.query_on_pos_strand() != b.query_on_pos_strand())
        && (a.query_on_pos_strand() == b.mate_on_pos_strand())
        && (a.mate_on_pos_strand() == b.query_on_pos_strand())
        && a.mapped_in_proper_pair()
        && b.mapped_in_proper_pair()
        && strcmp(amate_ref, bmate_ref) == 0;
}


//true if these entries are mate pairs and both unmapped
//used for storing a set of mate-paired reads when
bool AreUnmappedMatePairs(SamLine const& a, SamLine const& b)
{

    return
        AreSequencedMatePairs(a, b)
        && a.query_unmapped()
        && b.query_unmapped();
}



char const* strand_to_char(Strand strand)
{
    switch (strand)
    {
    case POS: return "POS"; break;
    case NEG: return "NEG"; break;
    case POS_NEG: return "POS_NEG"; break;
    }
    return "";
}



//3-way comparison operator for SAM by qname and flag
int SAM_cmp_qname_flag(SamLine const& a, SamLine const& b)
{

    return SAM_cmp_qname_flag_aux(a.qname, a.flag, b.qname, b.flag);
}


int SAM_cmp_qname_flag_aux(char const* qname1, int flag1,
                            char const* qname2, int flag2)
{
    int idcomp = strcmp(qname1, qname2);
    if (idcomp != 0)
    {
        return idcomp;
    }
    else
    {
        //subflag order
        int mask = SamFlags::PAIRED_IN_SEQUENCING 
            | SamFlags::FIRST_READ_IN_PAIR 
            | SamFlags::SECOND_READ_IN_PAIR;

        int a_subflag = flag1 & mask;
        int b_subflag = flag2 & mask;
        if (a_subflag < b_subflag)
        {
            return -1;
        }
        else if (a_subflag > b_subflag)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}


void print_sam_line(FILE * sam_fh,
                    char const* qname, int flag, 
                    char const* rname, size_t pos,
                    size_t mapq, char const* cigar,
                    char const* mrnm, size_t mpos,
                    int isize, char const* seq,
                    char const* qual,
                    char const* tag_string,
                    bool output_is_ones_based)
{
    int offset = output_is_ones_based ? 1 : 0;

    fprintf(sam_fh, "%s\t%i\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i\t%s\t%s",
            qname, flag, rname, pos + offset, mapq, cigar,
            mrnm, mpos + offset, isize, seq, qual);

    if (strlen(tag_string) > 0)
    {
        fprintf(sam_fh, "\t%s", tag_string);
    }
    fprintf (sam_fh, "\n");
}


//compute contig lengths from sam_fh.
//rewinds sam_fh to beginning
std::map<std::string, size_t> ContigLengths(FILE * sam_fh)
{
    SamLine * samline;
    std::map<std::string, size_t> contig_lengths;
    size_t total_length = 0;
    bool allow_absent_seq_qual = true;
    while ((samline = new SamLine(sam_fh, true, allow_absent_seq_qual)))
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
            contig_lengths[std::string(contig_name)] = contig_length;
            total_length += contig_length;
        }
        delete samline;
    }
    contig_lengths[std::string("*")] = 0;

    if (total_length == 0)
    {
        fprintf(stderr, "No SQ lines in header, or no SN / LN fields defining contigs\n");
        exit(1);
    }

    return contig_lengths;
}

CONTIG_OFFSETS ContigOffsets(std::map<std::string, size_t> const& contig_lengths)
{
    CONTIG_OFFSETS contig_offsets;

    size_t cumul_length = 0;
    for (std::map<std::string, size_t>::const_iterator cl = contig_lengths.begin();
         cl != contig_lengths.end(); ++cl)
    {
        char * contig_name = new char[(*cl).first.size() + 1];
        strcpy(contig_name, (*cl).first.c_str());
        contig_offsets[contig_name] = cumul_length;
        cumul_length += (*cl).second;
    }
    return contig_offsets;
}


//sets sam_fh to position of first data line
void SetToFirstDataLine(FILE ** sam_fh)
{
    SamLine * samline;
    rewind(*sam_fh);
    SAM_PARSE parse_flag;
    bool allow_absent_seq_qual = true;
    while ((samline = new SamLine(*sam_fh, true, allow_absent_seq_qual)))
    {
        if (samline->parse_flag == HEADER)
        {
            delete samline;
        }
        else
        {
            if (samline->line != NULL)
            {
                fseek(*sam_fh, - samline->bytes_in_line, std::ios::cur);
            }
            delete samline;
            break;
        }
    }
    return;
}


size_t flattened_position(char const* contig, size_t position, size_t flag,
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


//compute the flattened coordinate start position

//produce index consistent with ordering [rname, flag, pos]
size_t samline_position_align(char const* samline, CONTIG_OFFSETS const& contig_offsets)
{

    char rname[256];
    size_t pos;
    size_t flag;

    // qname flag rname pos mapq cigar mrnm mpos isize seq qual tags...
    int nfields_read = sscanf(samline, "%*s\t%i\t%s\t%zu", &flag, rname, &pos);
    
    if (nfields_read != 3)
    {
        fprintf(stderr, "samline_position_align: bad format for position-align sorting\n");
        exit(1);
    }

    CONTIG_OFFSETS::const_iterator dummy;
    return flattened_position(rname, pos, flag, contig_offsets, &dummy);
}


//compute the flattened coordinate start position.  only works for a read id that
//starts out numerically
size_t samline_read_id_flag(char const* samline, CONTIG_OFFSETS const& /* unused */)
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
size_t samline_position_min_align_guide(char const* samline, CONTIG_OFFSETS const& contig_offsets)
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

    CONTIG_OFFSETS::const_iterator guide_iter = contig_offsets.find(guide_chrom);

    if (guide_iter == contig_offsets.end())
    {
        fprintf(stderr, "This samline has its guide contig (%s) "
                "that is not in SAM contig index\n\n%s\n", 
                guide_chrom, samline);
        exit(1);
    }
    size_t guide_index = (*guide_iter).second + guide_pos;

    CONTIG_OFFSETS::const_iterator dummy;
    size_t align_index = 
        flattened_position(align_chrom, align_pos, flag, contig_offsets, &dummy);

    return std::min(guide_index, align_index);

}


char const* SamLine::mate_ref_name() const
{
    return strcmp(this->mrnm, "=") == 0 ? this->rname : this->mrnm;
}


void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh)
{
    SetToFirstDataLine(input_sam_fh);

    size_t header_length = ftell(*input_sam_fh);
    rewind(*input_sam_fh);

    char * header_buf = new char[header_length];
    fread(header_buf, 1, header_length, *input_sam_fh);
    fwrite(header_buf, 1, header_length, output_fh);
    delete header_buf;
    fflush(output_fh);
}

