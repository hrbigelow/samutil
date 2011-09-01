#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "sam_helper.h"
#include "sam_order.h"
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

//extern char const* AlignSpaceTag;

//serialize the data in a line.  string data will be packed in the
//'line' variable, with null's terminating each one.  non-string data
//will be stored separately

SAM_QNAME_FORMAT SamLine::sam_qname_format;
size_t (* SamLine::parse_fragment_id)(char const* qname) = &parse_fragment_id_zero;



void SamLine::SetGlobalFlags(SAM_QNAME_FORMAT _qname_format)
{
    SamLine::sam_qname_format = _qname_format;
    switch(SamLine::sam_qname_format)
    {
    case SAM_NUMERIC: SamLine::parse_fragment_id = &parse_fragment_id_numeric; break;
    case SAM_ILLUMINA: SamLine::parse_fragment_id = &parse_fragment_id_illumina; break;
    case SAM_CASAVA18: SamLine::parse_fragment_id = &parse_fragment_id_casava_1_8; break;
    case SAM_NON_INTERPRETED: SamLine::parse_fragment_id = &parse_fragment_id_zero; break;
    }
}


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
    isize(_isize), extra(NULL), extra_tag(NULL),
    flattened_pos(0), alignment_space(AlignSpaceMissing)
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

    this->fragment_id = SamLine::parse_fragment_id(this->qname);

    char as_tag[10];
    char as_type;
    if (this->has_tag(AlignSpaceTag, as_type, as_tag))
    {
        this->alignment_space = as_tag[0];
    }
}





SamLine::SamLine(FILE * seqfile, bool allow_absent_seq_qual) :
    bytes_in_line(0), line(NULL), qname(NULL), fragment_id(0), rname(NULL), 
    cigar(NULL), mrnm(NULL), seq(NULL), qual(NULL), tag_string(NULL), extra(NULL),
    extra_tag(NULL), flattened_pos(0), alignment_space(AlignSpaceMissing)
{
    char buf[4096 + 1];
    if (fgets(buf, 4096, seqfile) == NULL)
    {
        this->parse_flag = END_OF_FILE;
        return;
    }
    this->Init(buf, allow_absent_seq_qual);
}



SamLine::SamLine(char const* samline_string, bool allow_absent_seq_qual) :
    qname(NULL), fragment_id(0), rname(NULL), cigar(NULL), mrnm(NULL), 
    seq(NULL), qual(NULL), tag_string(NULL), extra(NULL), extra_tag(NULL),
    flattened_pos(0), alignment_space(AlignSpaceMissing)
{
    this->Init(samline_string, allow_absent_seq_qual);
}


//copy a SamLine according to the storage policy described in sam_helper.h
SamLine::SamLine(SamLine const& s) :
    bytes_in_line(s.bytes_in_line),
    parse_flag(s.parse_flag),
    fragment_id(s.fragment_id),
    flag(s.flag),
    pos(s.pos),
    mapq(s.mapq),
    mpos(s.mpos),
    isize(s.isize),
    flattened_pos(s.flattened_pos),
    alignment_space(s.alignment_space)
{
    this->line = new char[s.bytes_in_line];
    memcpy(this->line, s.line, s.bytes_in_line);
    this->qname = this->line + std::distance(s.line, s.qname);
    this->rname = this->line + std::distance(s.line, s.rname);
    this->mrnm = this->line + std::distance(s.line, s.mrnm);
    this->seq = this->line + std::distance(s.line, s.seq);
    this->qual = this->line + std::distance(s.line, s.qual);

    if (s.cigar == s.extra)
    {
        this->extra = new char[strlen(s.extra + 1)];
        strcpy(this->extra, s.extra);
        this->cigar = this->extra;
    }
    else
    {
        this->extra = NULL;
        this->cigar = this->line + std::distance(s.line, s.cigar);
    }

    if (s.tag_string == NULL && s.extra_tag == NULL)
    {
        this->tag_string = NULL;
        this->extra_tag = NULL;
    }
    else if (s.tag_string != NULL && s.extra_tag == NULL)
    {
        this->extra_tag = NULL;
        this->tag_string = this->line + std::distance(s.line, s.tag_string);
    }
    else
    {
        assert(s.tag_string == s.extra_tag);
        assert(s.extra_tag != NULL);

        this->extra_tag = new char[strlen(s.extra_tag + 1)];
        strcpy(this->extra_tag, s.extra_tag);
        this->tag_string = this->extra_tag;
    }
}


//samline_string must be terminated by newline
void SamLine::Init(char const* samline_string, bool allow_absent_seq_qual)
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

        this->fragment_id = SamLine::parse_fragment_id(this->qname);

        if (was_error)
        {
            fprintf(stderr, "SamLine::load_line: encountered badly formatted line.\n");
            fprintf(stderr, "Line:\n\n%s\n", this->line);
            this->parse_flag = PARSE_ERROR;
            return;
        }
        
        this->parse_flag = DATA_LINE;

        char as_tag[10];
        char as_type;
        if (this->has_tag(AlignSpaceTag, as_type, as_tag))
        {
            this->alignment_space = as_tag[0];
        }
    }
}


void SamLine::SetFlattenedPosition(CONTIG_OFFSETS const& contig_offsets,
                                   CONTIG_OFFSETS::const_iterator * contig_iter)
{
    //assume contig_iter is valid
    assert(*contig_iter != contig_offsets.end());
    //memoize the lookup, since we anticipate same lookups

    *contig_iter = strcmp((**contig_iter).first, this->rname) == 0
        ? *contig_iter 
        : contig_offsets.find(this->rname);
    
    if (*contig_iter == contig_offsets.end())
    {
        fprintf(stderr, "SetFlattenedPosition: error: rname %s (at %Zu) "
                "does not exist in provided contig index\n",
                this->rname, this->zero_based_pos());
        exit(1);
    }
    
    this->flattened_pos = (**contig_iter).second + this->zero_based_pos();
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

    int flag_to_print = flip_query_strand_flag ? 
        (this->flag ^ SamFlags::QUERY_ON_NEG_STRAND) : this->flag;

    size_t reported_pos = this->query_unmapped() ? 0 : this->ones_based_pos();
    size_t reported_mpos = this->mate_unmapped() ? 0 : this->mpos + 1;

    fprintf(seqfile, "%s\t%i\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i",
            this->qname, flag_to_print, this->rname, reported_pos,
            this->mapq, this->cigar, this->mrnm, reported_mpos,
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
    return parse_sam_tag(this->tag_string, tag, ":i:%i", default_if_missing, has_score);
}


/*
char SamLine::alignment_space(char default_if_missing, bool * has_space) const
{
    return parse_sam_tag(this->tag_string, AlignSpaceTag, ":%c", default_if_missing, has_space);
}
*/


bool SamLine::has_tag(char const* tag_code, char & type, char * value_string) const
{

    char const* tag_start;
    if (this->tag_string == NULL
        || (tag_start = strstr(this->tag_string, tag_code)) == NULL)
    {
        return false;
    }
    else
    {
        sscanf(tag_start, "%*[^:]:%c:%s", &type, value_string);
        return true;
    }
}


void SamLine::add_tag(char const* tag_code, char type, char const* value_string)
{

    char tag_and_value[20];
    sprintf(tag_and_value, "%s:%c:%s", tag_code, type, value_string);

    char tag_fragment[20];
    sprintf(tag_fragment, "%s:", tag_code);

    //1.  compute length of existing tag (or zero if not exists)
    char * tag_start = strstr(this->tag_string, tag_fragment);
    
    size_t old_tag_length = 0;
    if (tag_start != NULL)
    {
        char old_tag_and_value[20];
        sscanf(tag_start, "%s", old_tag_and_value);
        old_tag_length = strlen(old_tag_and_value);
    }
    size_t new_tag_length = strlen(tag_and_value);

    if (old_tag_length == new_tag_length)
    {
        //copy new value_string in place
        strncpy(tag_start, tag_and_value, new_tag_length);
    }
    else
    {
        //re-allocate and copy values
        char * new_extra_tag = new char[strlen(this->tag_string) + strlen(tag_and_value) + 2];
        
        if (tag_start != NULL)
        {
            //construct the new string in parts:
            size_t nbytes = std::distance(this->tag_string, tag_start);
            
            //1. everything before the start of the old tag (includes trailing space)
            strncpy(new_extra_tag, this->tag_string, nbytes);

            //2. the new tag (null-terminated)
            strcpy(new_extra_tag + nbytes, tag_and_value);

            //3. everything after the end of the old tag plus trailing tab.
            if (old_tag_length > 0)
            {
                strcat(new_extra_tag, tag_start + old_tag_length);
            }
        }
        else
        {
            strcpy(new_extra_tag, this->tag_string);
            strcat(new_extra_tag, "\t");
            strcat(new_extra_tag, tag_and_value);
        }

        if (this->extra_tag != NULL)
        {
            delete this->extra_tag;
        }
        this->extra_tag = new_extra_tag;
        this->tag_string = this->extra_tag;
    }
}


//true if samlines represent sequence mate pairs, regardless of
//whether they are mapped as such
bool AreSequencedMatePairs(SamLine const& a, SamLine const& b)
{
    return 
        a.fragment_id == b.fragment_id
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

    return
        AreSequencedMatePairs(a, b)
        && a.mpos == b.zero_based_pos()
        && b.mpos == a.zero_based_pos()
        && (a.query_on_pos_strand() != b.query_on_pos_strand())
        && (a.query_on_pos_strand() == b.mate_on_pos_strand())
        && (a.mate_on_pos_strand() == b.query_on_pos_strand())
        && a.mapped_in_proper_pair()
        && b.mapped_in_proper_pair()
        && strcmp(a.mate_ref_name(), b.mate_ref_name()) == 0;
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
/*
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
*/


/*
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
*/



//sets sam_fh to position of first data line
void SetToFirstDataLine(FILE ** sam_fh)
{
    SamLine * samline;
    rewind(*sam_fh);

    bool allow_absent_seq_qual = true;
    while ((samline = new SamLine(*sam_fh, allow_absent_seq_qual)))
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

