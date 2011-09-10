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
#include <algorithm>

// namespace SamFlags
// {
//     int const MULTI_FRAGMENT_TEMPLATE = 0x0001;
//     int const ALL_FRAGMENTS_MAPPED = 0x0002;
//     int const THIS_FRAGMENT_UNMAPPED = 0x0004;
//     int const NEXT_FRAGMENT_UNMAPPED = 0x0008;
//     int const THIS_FRAGMENT_ON_NEG_STRAND = 0x0010; //set when negative stranded
//     int const NEXT_FRAGMENT_ON_NEG_STRAND = 0x0020; //set when mate is negative stranded
//     int const FIRST_FRAGMENT_IN_TEMPLATE = 0x0040;
//     int const LAST_FRAGMENT_IN_TEMPLATE = 0x0080;
//     int const ALIGNMENT_NOT_PRIMARY = 0x0100;
//     int const FAILED_QUALITY_CHECK = 0x0200;
//     int const PCR_OR_OPTICAL_DUPLICATE = 0x0400;

//     //alternate rSAM definitions
//     int const TEMPLATE_ON_NEG_STRAND = 0x0010;

// };


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

char SamLine::expected_read_layout[256];

bool SamLine::expect_rsam_format;

size_t (* SamLine::parse_fragment_id)(char const* qname) = &parse_fragment_id_zero;



void SamLine::SetGlobalFlags(SAM_QNAME_FORMAT _qname_format,
                             char const* _expected_layout)
{
    SamLine::sam_qname_format = _qname_format;
    switch(SamLine::sam_qname_format)
    {
    case SAM_NUMERIC: SamLine::parse_fragment_id = &parse_fragment_id_numeric; break;
    case SAM_ILLUMINA: SamLine::parse_fragment_id = &parse_fragment_id_illumina; break;
    case SAM_CASAVA18: SamLine::parse_fragment_id = &parse_fragment_id_casava_1_8; break;
    case SAM_NON_INTERPRETED: SamLine::parse_fragment_id = &parse_fragment_id_zero; break;
    }

    strcpy(SamLine::expected_read_layout, _expected_layout);

    SamLine::expect_rsam_format = strcmp(SamLine::expected_read_layout, "") == 0;
}



SamLine::SamLine(SAM_PARSE _parse_flag,
                 char const* _qname, size_t _flagval, 
                 char const* _rname, size_t _pos,
                 size_t _mapq, char const* _cigar,
                 char const* _seq,
                 char const* _qual,
                 char const* _tag_string)
{
    char line[4096];
    if (_parse_flag == END_OF_FILE)
    {
        return;
    }
    SamFlag flag;
    flag.raw = _flagval;

    size_t reported_pos = flag.this_fragment_unmapped ? 0 : _pos + 1;

    if (SamLine::expect_rsam_format)
    {
        // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.SEQ.QUAL[.TAG[.TAG[.TAG...]]]

        sprintf(line, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%s",
                _qname, flag.raw, _rname, reported_pos,
                _mapq, _cigar, _seq, _qual);
    }
    else
    {
        // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
        char rnext[] = "=";
        size_t pnext = 0;
        int tlen = 0;
        sprintf(line, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i\t%s\t%s",
                _qname, flag.raw, _rname, reported_pos,
                _mapq, _cigar, rnext, pnext, tlen, _seq, _qual);
    }

    if (strlen(_tag_string) != 0)
    {
        strcat(line, "\t");
        strcat(line, _tag_string);
    }
    strcat(line, "\n");

    this->Init(line);
}


SamLine::SamLine(FILE * seqfile)
{
    char buf[4096 + 1];
    if (fgets(buf, 4096, seqfile) == NULL)
    {
        this->line = NULL;
        this->extra = NULL;
        this->extra_tag = NULL;
        this->parse_flag = END_OF_FILE;
        return;
    }
    this->Init(buf);
}



SamLine::SamLine(char const* samline_string)
{
    this->Init(samline_string);
}



SamLine::SamLine(SamLine const& s) :
    flattened_pos(s.flattened_pos)
{
    char * line = new char[s.bytes_in_line + 1];
    memcpy(line, s.line, s.bytes_in_line + 1);
    //replace all nulls with tabs
    char * l = line;
    while ((l = index(l, '\0')) != line + s.bytes_in_line - 1)
    {
        *l++ = '\t';
    }
    *l = '\n';

    this->Init(line);

    delete line;
}


//Merge two or more SAMlines that are part of the same fragment
SamLine::SamLine(SamLine const* samlines[], size_t n_lines, char const* read_layout)
    : parse_flag(samlines[0]->parse_flag),
      fragment_id(samlines[0]->fragment_id),
      pos(samlines[0]->pos),
      mapq(samlines[0]->mapq),
      pnext(0),
      tlen(0),
      flattened_pos(samlines[0]->flattened_pos),
      alignment_space(samlines[0]->alignment_space)
{
    flag.raw = 0;

    //assume at least one line
    assert(n_lines > 0);

    SamLine const* sl = samlines[0];

    //assume valid samlines, and that they are sorted by alignment position.
    //1 determine total space needed to store (bytes_in_line)

    // a. calculate sum of lengths of reads.
    size_t total_seq_length = strlen(sl->seq);

    // b. construct merged CIGAR
    char *merged_cigar = new char[1024];

    char cigar_op[16];
    Cigar::CIGAR_VEC prev_cigar = sl->flag.this_fragment_unmapped 
        ? Cigar::CIGAR_VEC({ Cigar::Unit(Cigar::Ops[Cigar::S], strlen(sl->seq))} )
        : Cigar::FromString(sl->cigar, 0);

    strcpy(merged_cigar, Cigar::ToString(prev_cigar).c_str());

    Cigar::CIGAR_VEC cur_cigar;

    for (size_t i = 1; i != n_lines; ++i)
    {
        assert(samlines[i]->alignment_space == this->alignment_space);

        total_seq_length += strlen(samlines[i]->seq);
        int64_t inter_seq_jump = 
            static_cast<int64_t>(samlines[i]->pos) 
            - (static_cast<int64_t>(samlines[i-1]->pos) 
               + static_cast<int64_t>(Cigar::Length(prev_cigar, true)));

        cur_cigar = samlines[i]->flag.this_fragment_unmapped 
            ? Cigar::CIGAR_VEC({ Cigar::Unit(Cigar::Ops[Cigar::S], strlen(samlines[i]->seq)) })
            : Cigar::FromString(samlines[i]->cigar, 0);
        
        sprintf(cigar_op, "%ziT", inter_seq_jump);
        strcat(merged_cigar, cigar_op);
        strcat(merged_cigar, Cigar::ToString(cur_cigar).c_str());

        prev_cigar = cur_cigar;
    }


    // c. merge tags appropriately (complicated!)  perhaps for now, we
    // simply can catenate tag strings together in fragment alignment
    // order. then, when printing, de-duplicate them...

    // for now, just take the tag string of the first record as the
    // string for the whole one
    
    
    // d. with all data, calculate bytes_in_line and allocate
    this->bytes_in_line = 
        strlen(sl->qname)
        + strlen(sl->rname)
        + strlen(merged_cigar)
        + (2 * total_seq_length)
        + strlen(sl->tag_string)
        + 7; // space for null terminators

    this->line = new char[this->bytes_in_line + 1];

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

    this->tag_string = this->qual + strlen(this->qual) + 1;
    strcpy(this->tag_string, sl->tag_string);

    this->extra = NULL;
    this->extra_tag = NULL;

    this->flag.raw = 0;

    this->flag.is_rsam_format = 1;
    this->flag.num_fragments_in_template = n_lines;
    for (size_t i = 0; i != n_lines; ++i)
    {
        this->flag.read_layout |= (read_layout[i] == 'f') ? (1<<i) : 0;
    }

    //these flags are dominant.  If they are set in any one of the
    //reads, they should be set in the merged record
    for (size_t i = 0; i != n_lines; ++i)
    {
        this->flag.failed_quality_check |= samlines[i]->flag.failed_quality_check;
        this->flag.pcr_or_optical_duplicate |= samlines[i]->flag.pcr_or_optical_duplicate;
    }    

    // c. initialize flag. assume (but do not check) that these flags
    // have the same values for all records on this template.
    this->flag.all_fragments_mapped |= sl->flag.all_fragments_mapped;
    this->flag.alignment_not_primary |= sl->flag.alignment_not_primary;

    // reinterpret 0x10 as 'fragment orientation'
    // forward fragments will have left-most read be 'first-in-template'
    this->flag.this_fragment_on_neg_strand |= samlines[0]->flag.first_fragment_in_template;
    
    delete merged_cigar;
}


//samline_string must be terminated by newline
void SamLine::Init(char const* samline_string)
{

    this->extra = NULL;
    this->extra_tag = NULL;

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

        int off[5];
        int tab[5]; // mark positions just before tabs

        char * rest_of_line;

        if (SamLine::expect_rsam_format)
        {
            //QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
            int num_fields =
                sscanf(this->line, 
                       "%n%*s%n\t" "%zu\t" "%n%*s%n\t" "%zu\t" "%zu\t"
                       "%n%*s%n\t" "%n",
                       off, tab, &this->flag.raw, off + 1, tab + 1, &this->pos, &this->mapq, 
                       off + 2, tab + 2, off + 3);

            if (num_fields != 3
                || ! std::is_sorted(off, off + 4)
                || ! std::is_sorted(tab, tab + 3))
            {
                was_error = true;
            }
            else
            {
                this->qname = this->line + off[0];
                this->rname = this->line + off[1];
                this->cigar = this->line + off[2];
                rest_of_line = this->line + off[3];

                this->rnext = NULL;
                this->pnext = 0;
                
                this->line[tab[0]] = '\0';
                this->line[tab[1]] = '\0';
                this->line[tab[2]] = '\0';
            }
        }
        else
        {
            //QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
            int num_fields =
                sscanf(this->line, 
                       "%n%*s%n\t" "%zu\t" "%n%*s%n\t" "%zu\t" "%zu\t" 
                       "%n%*s%n\t" "%n%*s%n\t" 
                       "%zu\t" "%i\t" "%n",
                       off, tab, &this->flag.raw, off + 1, tab + 1, &this->pos, &this->mapq, 
                       off + 2, tab + 2, off + 3, tab + 3,
                       &this->pnext, &this->tlen, off + 4);

            if (num_fields != 5
                || ! std::is_sorted(off, off + 5)
                || ! std::is_sorted(tab, tab + 4))
            {
                was_error = true;
            }
            else
            {
                this->qname = this->line + off[0];
                this->rname = this->line + off[1];
                this->cigar = this->line + off[2];
                this->rnext = this->line + off[3];
                rest_of_line = this->line + off[4];
                
                this->line[tab[0]] = '\0';
                this->line[tab[1]] = '\0';
                this->line[tab[2]] = '\0';
                this->line[tab[3]] = '\0';
            }
        }

        //rest_of_line must be one of:
        //SEQ.QUAL.TAGS
        //SEQ.QUAL
        //.TAGS
        //.

        int rest_off[3];

        if (! was_error)
        {
            //handles cases SEQ.QUAL.TAGS and SEQ.QUAL
            if (sscanf(rest_of_line, "%n%*s%n\t" "%n%*s", rest_off, tab, rest_off + 1) == 0)
            {
                this->seq = rest_of_line + rest_off[0];
                this->qual = rest_of_line + rest_off[1];
                rest_of_line[tab[0]] = '\0'; // truncates 'seq'
                //search for a 'tab' in qual field.
                char *qual_end = strchr(this->qual, '\t');
                if (qual_end == NULL)
                {
                    this->tag_string = NULL;
                }
                else
                {
                    *qual_end = '\0';
                    this->tag_string = qual_end + 1;
                }
            }
            
            //handles cases . and .TAGS
            else if (sscanf(rest_of_line, "\t%n", rest_off) == 0)
            {
                //seq and qual are zero length, there are no tags
                this->seq = rest_of_line;
                this->qual = rest_of_line + 1;
                this->tag_string = rest_of_line + 2;
                this->seq[0] = '\0';
                this->qual[0] = '\0';
            }
            
            else
            {
                was_error = true;
            }
        }

        if (was_error)
        {
            char const* format = SamLine::expect_rsam_format ? "rSAM" : "SAM";
            fprintf(stderr, "Error: SamLine: Bad %s format:\n\n%s\n",
                    format, this->line);
            this->parse_flag = PARSE_ERROR;
            return;
        }
        
        if (! this->flag.this_fragment_unmapped)
        {
            --this->pos;
        }
        if (! this->flag.is_rsam_format && ! this->flag.next_fragment_unmapped)
        {
            --this->pnext;
        }            

        this->fragment_id = SamLine::parse_fragment_id(this->qname);

        this->parse_flag = DATA_LINE;

        this->flattened_pos = 0;

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

    else if (! this->flag.is_rsam_format)
    {
        size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();

        fprintf(seqfile, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i",
                this->qname, flag.raw, this->rname, reported_pos,
                this->mapq, this->cigar, this->rnext, this->pnext, this->tlen);

        if (this->tag_string != NULL)
        {
            fprintf(seqfile, "\t%s", this->tag_string);
        }
        fprintf (seqfile, "\n");

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

    SamFlag flag = this->flag;

    int layout_index = 0;
    int this_fragment_layout;
    int next_fragment_layout;

    size_t num_fragments = flag.num_fragments_in_template;

    size_t read_index = 0;

    while (start != cigar_v.end())
    {
        this_fragment_layout = this->flag.read_layout & (1<<layout_index);
        next_fragment_layout = this->flag.read_layout & (1<<((layout_index + 1) % num_fragments));

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
        }
        else
        {
            //multi-fragment template
            int64_t t_oplen = (*end).length;
            pnext = cur_read_start + read_length + t_oplen;
            ++end;
            flag.multi_fragment_template = 1;
            
            //fragment_layout is fragment-to-template.
            //template_layout is template-to-reference.
            //we want to calculate fragment-to-reference.
            flag.this_fragment_on_neg_strand =
                (this_fragment_layout == this->flag.template_layout) 
                ? LAYOUT_FORWARD
                : LAYOUT_REVERSE;

            flag.next_fragment_on_neg_strand =
                (next_fragment_layout == this->flag.template_layout) 
                ? LAYOUT_FORWARD
                : LAYOUT_REVERSE;
        }
        start = end;

        int64_t used_tlen = read_index == 0 
            ? tlen 
            : (read_index == num_fragments - 1 ? -tlen : 0);

        fprintf(seqfile, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%zi",
                this->qname, flag.raw, this->rname, cur_read_start,
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

        ++layout_index;
    }

}



void SamLine::print_rsam(FILE * seqfile) const
{

    if (this->parse_flag == HEADER)
    {
        fprintf(seqfile, "%s\n", this->line);
        return;
    }

    size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();

    fprintf(seqfile, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%s",
            this->qname, this->flag.raw, this->rname, reported_pos,
            this->mapq, this->cigar, this->seq, this->qual);

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
        && (a.flag.first_fragment_in_template != b.flag.first_fragment_in_template)
        && a.flag.multi_fragment_template
        && b.flag.multi_fragment_template;
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
        && (a.flag.this_fragment_on_neg_strand == b.flag.next_fragment_on_neg_strand)
        && (a.flag.next_fragment_on_neg_strand == b.flag.this_fragment_on_neg_strand)
        && a.flag.all_fragments_mapped
        && b.flag.all_fragments_mapped
        && strcmp(a.next_fragment_ref_name(), b.next_fragment_ref_name()) == 0;
}


//true if these entries are mate pairs and both unmapped
//used for storing a set of mate-paired reads when
bool AreUnmappedMatePairs(SamLine const& a, SamLine const& b)
{

    return
        AreSequencedMatePairs(a, b)
        && a.flag.this_fragment_unmapped
        && b.flag.this_fragment_unmapped;
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
    rewind(*sam_fh);

    //set just after first newline that doesn't begin with '@'
    char c;
    size_t dummy;
    char *line = NULL;

    while ((c = fgetc(*sam_fh)) == '@')
    {
        ungetc(c, *sam_fh);
        getline(&line, &dummy, *sam_fh);
    }
    ungetc(c, *sam_fh);

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



