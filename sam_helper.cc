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
#include <functional>
#include <string>

namespace SamFlags
{
    int const MULTI_FRAGMENT_TEMPLATE = 0x0001;
    int const ALL_FRAGMENTS_MAPPED = 0x0002;
    int const THIS_FRAGMENT_UNMAPPED = 0x0004;
    int const NEXT_FRAGMENT_UNMAPPED = 0x0008;
    int const THIS_FRAGMENT_ON_NEG_STRAND = 0x0010; //set when negative stranded
    int const NEXT_FRAGMENT_ON_NEG_STRAND = 0x0020; //set when mate is negative stranded
    int const FIRST_FRAGMENT_IN_TEMPLATE = 0x0040;
    int const LAST_FRAGMENT_IN_TEMPLATE = 0x0080;
    int const ALIGNMENT_NOT_PRIMARY = 0x0100;
    int const FAILED_QUALITY_CHECK = 0x0200;
    int const PCR_OR_OPTICAL_DUPLICATE = 0x0400;

    //alternate rSAM definitions
    int const TEMPLATE_ON_NEG_STRAND = 0x0010;

};


bool eqstr::operator()(const char* s1, const char* s2) const
{
    return strcmp(s1, s2) == 0;
}


//a hashing function.  Expects that k will either start or end with a string of digits.
//if it does, the integer value of those digits is used.
//Otherwise, if it is '*', zero is used.,
//Otherwise, the default hash of its string value is used.
size_t to_integer::operator()(char const* k) const
{
    size_t val = 0;
    char dummy; //used as a test whether sscanf reached the end of the string
    size_t dummy_int;

    if (sscanf(k, "%*[^0-9]%zu%c", &val, &dummy) == 1 
        || sscanf(k, "%zu%*[^0-9]%zu", &val, &dummy_int) == 1)
    {
        return val;
    }
    else if (k[0] == '*')
    {
        return 0;
    }
    else
    {
        return static_cast<STRING_HASH>(*this)(std::string(k, strlen(k)));
    }
}

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
                 char const* _rnext, size_t _pnext,
                 int _isize, char const* _seq,
                 char const* _qual,
                 char const* _tag_string) :
    parse_flag(_parse_flag), flag(_flag), pos(_pos), 
    mapq(_mapq), pnext(_pnext), 
    isize(_isize), extra(NULL), extra_tag(NULL),
    flattened_pos(0), alignment_space(AlignSpaceMissing)
{
    
    assert(strlen(_qual) < 10000);
    assert(strlen(_seq) < 10000);

    char const* src[] = { 
        _qname, _rname, _cigar, _rnext, _seq, _qual, _tag_string
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
        &this->rnext, &this->seq, &this->qual,
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
    cigar(NULL), rnext(NULL), seq(NULL), qual(NULL), tag_string(NULL), extra(NULL),
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
    qname(NULL), fragment_id(0), rname(NULL), cigar(NULL), rnext(NULL), 
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
    pnext(s.pnext),
    isize(s.isize),
    flattened_pos(s.flattened_pos),
    alignment_space(s.alignment_space)
{
    this->line = new char[s.bytes_in_line];
    memcpy(this->line, s.line, s.bytes_in_line);
    this->qname = this->line + std::distance(s.line, s.qname);
    this->rname = this->line + std::distance(s.line, s.rname);
    this->rnext = this->line + std::distance(s.line, s.rnext);
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



//Merge two or more SAMlines that are part of the same fragment
SamLine::SamLine(SamLine const* samlines[], size_t n_lines, char const* read_layout)
    : parse_flag(samlines[0]->parse_flag),
      fragment_id(samlines[0]->fragment_id),
      flag(0),
      pos(samlines[0]->pos),
      mapq(samlines[0]->mapq),
      pnext(0),
      isize(0),
      flattened_pos(samlines[0]->flattened_pos)
{
    //assume at least one line
    assert(n_lines > 0);

    SamLine const* sl = samlines[0];

    //assume valid samlines, and that they are sorted by alignment position.
    //1 determine total space needed to store (bytes_in_line)

    // a. calculate sum of lengths of reads.
    size_t total_seq_length = strlen(sl->seq);

    // b. construct CIGAR with excess memory.  (simply catenate together with intervening 'H')
    char *merged_cigar = new char[1024];
    strcpy(merged_cigar, sl->cigar);

    char cigar_op[16];
    Cigar::CIGAR_VEC prev_cigar = Cigar::FromString(sl->cigar, 0);

    for (size_t i = 1; i != n_lines; ++i)
    {
        total_seq_length += strlen(samlines[i]->seq);
        int64_t inter_seq_jump = 
            static_cast<int64_t>(samlines[i]->pos) 
            - (static_cast<int64_t>(samlines[i-1]->pos) 
               + static_cast<int64_t>(Cigar::Length(prev_cigar, true)));

        sprintf(cigar_op, "%ziH", inter_seq_jump);
        strcat(merged_cigar, cigar_op);
        strcat(merged_cigar, samlines[i]->cigar);
    }


    // c. merge tags appropriately (complicated!)  perhaps for now, we
    // simply can catenate tag strings together in fragment alignment
    // order. then, when printing, de-duplicate them...
    
    // d. with all data, calculate bytes_in_line and allocate
    this->bytes_in_line = 
        strlen(sl->qname)
        + strlen(sl->rname)
        + strlen(merged_cigar)
        + (2 * total_seq_length)
        + strlen(read_layout)
        + 7; // space for null terminators

    this->line = new char[this->bytes_in_line];

    //2 initialize all remaining fields and pointers

    // a. copy qname, rname, cigar, all seqs, all quals, read layout, and tags to 'line'
    this->qname = this->line;
    strcpy(this->qname, sl->qname);

    this->rname = this->qname + strlen(this->qname) + 1;
    strcpy(this->rname, sl->rname);

    this->cigar = this->rname + strlen(this->rname) + 1;
    strcpy(this->cigar, merged_cigar);

    this->seq = this->cigar + strlen(this->cigar) + 1;
    this->seq[0] = '\0';
    for (size_t i = 0; i != n_lines; ++i)
    {
        strcat(this->seq, samlines[i]->seq);
    }    

    this->qual = this->seq + strlen(this->seq) + 1;
    this->qual[0] = '\0';
    for (size_t i = 0; i != n_lines; ++i)
    {
        strcat(this->qual, samlines[i]->qual);
    }    

    this->read_layout = this->qual + strlen(this->qual) + 1;
    strcpy(this->read_layout, read_layout);

    this->flag = 0;
    //these flags are dominant.  If they are set in any one of the
    //reads, they should be set in the merged record
    for (size_t i = 0; i != n_lines; ++i)
    {
        this->flag |= (samlines[i]->flag & SamFlags::FAILED_QUALITY_CHECK);
        this->flag |= (samlines[i]->flag & SamFlags::PCR_OR_OPTICAL_DUPLICATE);
    }    

    // c. initialize flag. assume (but do not check) that these flags
    // have the same values for all records on this template.
    this->flag |=
        (sl->flag & SamFlags::ALL_FRAGMENTS_MAPPED)
        | (sl->flag & SamFlags::ALIGNMENT_NOT_PRIMARY);

    // reinterpret 0x10 as 'fragment orientation'
    // forward fragments will have left-most read be 'first-in-template'
    this->flag |= (samlines[0]->first_fragment_in_template() ? 0 : SamFlags::THIS_FRAGMENT_ON_NEG_STRAND);
    
    delete merged_cigar;
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
                   &this->flag, &this->pos, &this->mapq, &this->pnext, &this->isize);

        if (num_numeric_fields != 5)
        {
            was_error = true;
        }

        //now, replace tabs with nulls, and set addresses
        char * ptr = this->line;

        char * dummy = NULL;

        char ** trg[] = { 
            &this->qname, &dummy, &this->rname, &dummy, &dummy, &this->cigar, &this->rnext,
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

        if (*trg_ptr != &this->tag_string && ! (allow_absent_seq_qual && seq_absent))
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
        

        if (! this->this_fragment_unmapped())
        {
            --this->pos;
        }
        if (! this->next_fragment_unmapped())
        {
            --this->pnext;
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


// 
void SamLine::print_sam(FILE * seqfile) const
{

    if (this->parse_flag == HEADER)
    {
        fprintf(seqfile, "%s\n", this->line);
        return;
    }

    size_t cur_seq_start = 0;
    int64_t cur_read_start = this->pos;
    int64_t pnext;
    char rnext[2] = "=";
    
    Cigar::CIGAR_VEC cigar_v = Cigar::FromString(this->cigar, 0);

    Cigar::CIGAR_ITER start = cigar_v.begin();
    Cigar::CIGAR_ITER end = cigar_v.begin();

    int64_t tlen = Cigar::Length(cigar_v, false);

    char cigar_substring[1024];

    int flag = this->flag;

    char const* this_read_layout = this->read_layout;
    char const* next_read_layout = this->read_layout;

    size_t num_reads = strlen(this->read_layout);
    size_t read_index = 0;

    while (start != cigar_v.end())
    {
        while (end != cigar_v.end() 
               && (*end).op.code != Cigar::T
               && (*end).op.code != Cigar::U)
        {
            ++end;
        }
        //now 'end' points to the real end or to a 'T' operator
        size_t read_length = Cigar::Length(start, end, false);

        Cigar::ToString(start, end, cigar_substring);

        if (end == cigar_v.end())
        {
            //we're at the end of the series of reads.
            //set pnext to point back to the first pos
            pnext = this->pos;
            next_read_layout = this->read_layout;
        }
        else
        {
            //we've found a 'T' or 'U' operator.  So, this is a multiple-fragment template.
            int64_t t_oplen = (*end).length;
            pnext = cur_read_start + read_length + t_oplen;
            ++end;
            flag ^= SamFlags::MULTI_FRAGMENT_TEMPLATE;
            
            if ((*this_read_layout == 'f') == this->template_on_pos_strand()) 
            {
                //
                flag &= ~(SamFlags::THIS_FRAGMENT_ON_NEG_STRAND);
            }
            else
            {
                flag |= SamFlags::THIS_FRAGMENT_ON_NEG_STRAND;
            }
            if ((*next_read_layout == 'f') == this->template_on_pos_strand())
            {
                flag &= ~(SamFlags::NEXT_FRAGMENT_ON_NEG_STRAND);
            }
            else
            {
                flag |= SamFlags::NEXT_FRAGMENT_ON_NEG_STRAND;
            }
        }
        start = end;

        int64_t used_tlen = read_index == 0 
            ? tlen 
            : (read_index == num_reads - 1 ? -tlen : 0);

        fprintf(seqfile, "%s\t%i\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%zi",
                this->qname, flag, this->rname, cur_read_start,
                this->mapq, cigar_substring, rnext, pnext, used_tlen);
        
        if (this->seq != NULL)
        {
            assert(this->qual != NULL);
            fprintf(seqfile, "\t");
            fwrite(this->seq + cur_seq_start, 1, read_length, seqfile);
            fprintf(seqfile, "\t");
            fwrite(this->qual + cur_seq_start, 1, read_length, seqfile);
            cur_seq_start += read_length;
        }

        //now about the tags.
        if (this->tag_string != NULL)
        {
            fprintf(seqfile, "\t%s", this->tag_string);
        }
        fprintf (seqfile, "\n");

        ++this_read_layout;
        ++next_read_layout;
    }

}



void SamLine::print_rsam(FILE * seqfile) const
{

    if (this->parse_flag == HEADER)
    {
        fprintf(seqfile, "%s\n", this->line);
        return;
    }

    size_t reported_pos = this->this_fragment_unmapped() ? 0 : this->ones_based_pos();

    fprintf(seqfile, "%s\t%i\t%s\t%Zu\t%Zu\t%s\t%s\t%s\t%s",
            this->qname, this->flag, this->rname, reported_pos,
            this->mapq, this->cigar, this->read_layout,
            this->seq, this->qual);

    if (this->tag_string != NULL)
    {
        fprintf(seqfile, "\t%s", this->tag_string);
    }
    fprintf (seqfile, "\n");
}


void SamLine::print(FILE * seqfile, bool print_rsam) const
{
    if (print_rsam)
    {
        this->print_rsam(seqfile);
    }
    else
    {
        this->print_sam(seqfile);
    }
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
        && (a.first_fragment_in_template() != b.first_fragment_in_template())
        && a.multi_fragment_template()
        && b.multi_fragment_template();
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
        && a.pnext == b.zero_based_pos()
        && b.pnext == a.zero_based_pos()
        && (a.this_fragment_on_pos_strand() != b.this_fragment_on_pos_strand())
        && (a.this_fragment_on_pos_strand() == b.next_fragment_on_pos_strand())
        && (a.next_fragment_on_pos_strand() == b.this_fragment_on_pos_strand())
        && a.all_fragments_mapped()
        && b.all_fragments_mapped()
        && strcmp(a.next_fragment_ref_name(), b.next_fragment_ref_name()) == 0;
}


//true if these entries are mate pairs and both unmapped
//used for storing a set of mate-paired reads when
bool AreUnmappedMatePairs(SamLine const& a, SamLine const& b)
{

    return
        AreSequencedMatePairs(a, b)
        && a.this_fragment_unmapped()
        && b.this_fragment_unmapped();
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


char const* SamLine::next_fragment_ref_name() const
{
    return strcmp(this->rnext, "=") == 0 ? this->rname : this->rnext;
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



