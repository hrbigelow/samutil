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


SamFilter::SamFilter(char const* _tag_filter,
                     size_t _min_mapq,
                     size_t _max_stratum_rank,
                     size_t _max_stratum_size,
                     char const* _align_space_filt) :

    min_mapping_quality(_min_mapq),
    max_stratum_rank(_max_stratum_rank),
    max_stratum_size(_max_stratum_size),
    alignment_space_filter(_align_space_filt)
{ 
    // set masking flags
    this->comb.set_raw(0);
    this->values.set_raw(0);
    
    for (char const* c = _tag_filter; *c != '\0'; ++c)
    {
        bool v = isupper(*c);
        
        switch (toupper(*c))
        {
        case 'M': this->comb.all_fragments_mapped = true; this->values.all_fragments_mapped = v; break;
        case 'P': this->comb.alignment_not_primary = true; this->values.alignment_not_primary = ! v; break;
        case 'Q': this->comb.failed_quality_check = true; this->values.failed_quality_check = ! v; break;
        case 'D': this->comb.pcr_or_optical_duplicate = true; this->values.pcr_or_optical_duplicate = ! v; break;
        default: 
            fprintf(stderr, "Error: SamFilter: tag filter '%s' has bad content. Should be\n"
                    "[Mm][Pp][Qq][Dd]\n", _tag_filter);
            exit(1);
            break;
        }
    }
}



bool SamFilter::pass(SamLine const* samline) const
{
    return (samline->flag.get_raw() & this->comb.get_raw()) == this->values.get_raw()
        
        && samline->mapq >= this->min_mapping_quality
        
        && ((! samline->tags.stratum_rank_present)
            || samline->tags.stratum_rank <= this->max_stratum_rank)

        && ((! samline->tags.stratum_size_present)
            || samline->tags.stratum_size <= this->max_stratum_size)

        && ((! samline->tags.alignment_space_present)
            || this->alignment_space_filter == NULL
            || strchr(this->alignment_space_filter, samline->tags.alignment_space) != NULL);

}


size_t SamFlag::get_raw() const
{
    size_t raw = static_cast<size_t>
        (this->multi_fragment_template % 2
         | (this->all_fragments_mapped % 2)<<1
         | (this->this_fragment_unmapped % 2)<<2
         | (this->next_fragment_unmapped % 2)<<3
         | (this->this_fragment_on_neg_strand % 2)<<4
         | (this->next_fragment_on_neg_strand % 2)<<5
         | (this->first_fragment_in_template % 2)<<6
         | (this->last_fragment_in_template % 2)<<7
         | (this->alignment_not_primary % 2)<<8
         | (this->failed_quality_check % 2)<<9
         | (this->pcr_or_optical_duplicate % 2)<<10
         | (this->template_layout % 2)<<11
         )
        | static_cast<size_t>(this->is_rsam_format % (1<<8))<<16
        | static_cast<size_t>(this->num_fragments_in_template % (1<<8))<<24
        | static_cast<size_t>(this->read_layout % (1L<<32))<<32;

    return raw;
}


void SamFlag::set_raw(size_t raw)
{
    this->multi_fragment_template = static_cast<unsigned int>(raw % 2);
    this->all_fragments_mapped = static_cast<unsigned int>((raw>>1) % 2);
    this->this_fragment_unmapped = static_cast<unsigned int>((raw>>2) % 2);
    this->next_fragment_unmapped = static_cast<unsigned int>((raw>>3) % 2);
    this->this_fragment_on_neg_strand = static_cast<unsigned int>((raw>>4) % 2);
    this->next_fragment_on_neg_strand = static_cast<unsigned int>((raw>>5) % 2);
    this->first_fragment_in_template = static_cast<unsigned int>((raw>>6) % 2);
    this->last_fragment_in_template = static_cast<unsigned int>((raw>>7) % 2);
    this->alignment_not_primary = static_cast<unsigned int>((raw>>8) % 2);
    this->failed_quality_check = static_cast<unsigned int>((raw>>9) % 2);
    this->pcr_or_optical_duplicate = static_cast<unsigned int>((raw>>10) % 2);
    this->template_layout = static_cast<unsigned int>((raw>>11) % 2);

    this->is_rsam_format = static_cast<unsigned int>((raw>>16) % (1<<8));
    this->num_fragments_in_template = static_cast<unsigned int>((raw>>24) % (1<<8));
    this->read_layout = static_cast<unsigned int>((raw>>32) % (1<<8));

}


size_t SamTag::get_raw() const
{
    size_t raw =
        (this->raw_score % (1<<12))
        | static_cast<size_t>(this->stratum_rank % (1<<10))<<12
        | static_cast<size_t>(this->stratum_size % (1<<12))<<22
        | static_cast<size_t>(this->alignment_space)<<34
        | static_cast<size_t>(this->raw_score_present % (1<<1))<<42
        | static_cast<size_t>(this->alignment_space_present % (1<<1))<<43
        | static_cast<size_t>(this->stratum_rank_present % (1<<1))<<44
        | static_cast<size_t>(this->stratum_size_present % (1<<1))<<45
        | static_cast<size_t>(this->raw_score_tag[0])<<46
        | static_cast<size_t>(this->raw_score_tag[1])<<54;

    return raw;
}


void SamTag::set_raw(size_t raw)
{
    this->raw_score = static_cast<uint16_t>(raw % (1<<12));
    this->stratum_rank = static_cast<uint16_t>((raw>>12) % (1<<10));
    this->stratum_size = static_cast<uint16_t>((raw>>22) % (1<<12));
    this->alignment_space = static_cast<char>((raw>>34) % (1<<8));
    this->raw_score_present = static_cast<bool>((raw>>42) % (1<<1));
    this->alignment_space_present = static_cast<bool>((raw>>43) % (1<<1));
    this->stratum_rank_present = static_cast<bool>((raw>>44) % (1<<1));
    this->stratum_size_present = static_cast<bool>((raw>>45) % (1<<1));

    this->raw_score_tag[0] = static_cast<char>((raw>>46) % (1<<8));
    this->raw_score_tag[1] = static_cast<char>((raw>>54) % (1<<8));
}

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


char base_to_complement[] =
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XTXGXXXCXXXXXXNXXXXXAXXXXXXXXXXX"
    "XtXgXXXcXXXXXXnXXXXXaXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";


// copies the reverse complement of a range defined by [begin, end)
// into rcomp.  rcomp must not overlap [begin, end) range
void reverse_comp(char const* begin, char const* end, char * rcomp)
{

    while (begin != end)
    {
        --end;
        *rcomp = base_to_complement[static_cast<int>(*end)];
        ++rcomp;
    }
}


// outputs the reverse complement of a range defined by [begin, end)
// into rcomp.  rcomp may be equal to begin or not, but it must contain
// enough space
void reverse_comp_inplace(char const* begin, char const* end)
{
    
    char * comp;
    for (comp = begin; comp != end; ++comp)
    {
        *comp = base_to_complement[static_cast<int>(*comp)];
    }
    std::reverse(begin, end);
}



//
void init_tags_from_string(char const* tag_string, 
                           char const* raw_score_tag,
                           size_t worst_raw_score, 
                           SamTag * tags,
                           char ** extra_tags)
{

    // This function should only be used when initializing from SAM
    // input, not rSAM input.
    assert(! SamLine::expect_rsam_format);

    char tmp_buf[1024];
    char * tmp_buf_ptr = tmp_buf;
    tmp_buf_ptr[0] = '\0';

    char const* tag = tag_string;
    char tag_name[3];
    char tag_type;
    char tag_value[1024];
    int advance;

    // set default values.  
    (*tags).raw_score = 0;
    (*tags).stratum_rank = 0;
    (*tags).stratum_size = 0;
    (*tags).alignment_space = '\0';
    (*tags).raw_score_present = false;
    (*tags).alignment_space_present = false;
    (*tags).stratum_rank = false;
    (*tags).stratum_size_present = false;
    strcpy((*tags).raw_score_tag, "--");
    
    
    while (sscanf(tag, " %[^:]:%c:%s%n", tag_name, &tag_type, tag_value, &advance) == 3)
    {
        if (strcmp(tag_name, raw_score_tag) == 0)
        {
            strcpy((*tags).raw_score_tag, raw_score_tag);
            (*tags).raw_score = atoi(tag_value);
            (*tags).raw_score_present = true;
        }
        else if (strcmp(tag_name, AlignSpaceTag) == 0)
        {
            (*tags).alignment_space = tag_value[0];
            (*tags).alignment_space_present = true;
        }
        else if (strcmp(tag_name, StratumRankTag) == 0)
        {
            (*tags).stratum_rank = atoi(tag_value);
            (*tags).stratum_rank_present = true;
        }
        else if (strcmp(tag_name, StratumSizeTag) == 0)
        {
            (*tags).stratum_size = atoi(tag_value);
            (*tags).stratum_size_present = true;
        }
        else
        {
            if (tmp_buf_ptr != tmp_buf)
            {
                //we've added at least one tag.  tab-separate any additional ones
                tmp_buf_ptr += sprintf(tmp_buf_ptr, "\t");
            }
            tmp_buf_ptr += sprintf(tmp_buf_ptr, "%s:%c:%s", tag_name, tag_type, tag_value);
        }
        tag += advance;
    }

    if (strlen(tmp_buf) > 0)
    {
        *extra_tags = new char[strlen(tmp_buf) + 1];
        strcpy(*extra_tags, tmp_buf);
    }
    else
    {
        *extra_tags = NULL;
    }

}



//extern char const* AlignSpaceTag;

//serialize the data in a line.  string data will be packed in the
//'line' variable, with null's terminating each one.  non-string data
//will be stored separately

SAM_QNAME_FORMAT SamLine::sam_qname_format;

char SamLine::expected_read_layout[256];
bool SamLine::expect_rsam_format;
char SamLine::raw_score_tag[3];
size_t SamLine::worst_fragment_score;
bool SamLine::initialized = false;
bool SamLine::retain_qname_string = false;

size_t (* SamLine::parse_fragment_id)(char const* qname) = &parse_fragment_id_zero;



void SamLine::SetGlobalFlags(SAM_QNAME_FORMAT _qname_format,
                             char const* _expected_layout,
                             char const* _raw_score_tag,
                             size_t _worst_fragment_score,
                             bool _retain_qname_string)
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
    strcpy(SamLine::raw_score_tag, _raw_score_tag);
    SamLine::worst_fragment_score = _worst_fragment_score;

    SamLine::expect_rsam_format = strcmp(SamLine::expected_read_layout, "") == 0;
    SamLine::retain_qname_string = _retain_qname_string;
    SamLine::initialized = true;
}



SamLine::SamLine(SAM_PARSE _parse_flag,
                 char const* _qname, size_t _flagval, 
                 char const* _rname, size_t _pos,
                 size_t _mapq, char const* _cigar,
                 size_t _tagval,
                 char const* _tag_string)
{
    char line[4096];
    if (_parse_flag == END_OF_FILE)
    {
        return;
    }
    SamFlag flag;
    flag.set_raw(_flagval);

    size_t reported_pos = flag.this_fragment_unmapped ? 0 : _pos + 1;

    if (SamLine::expect_rsam_format)
    {
        // FID.FLAG.RNAME.POS.MAPQ.CIGAR.TAGRAW[.TAG[.TAG[.TAG...]]]
        size_t fragment_id = this->parse_fragment_id(_qname);
        if (fragment_id == QNAME_FORMAT_ERROR)
        {
            fprintf(stderr, "Error: samline_fragment_id: bad qname format in this qname\n%s\n\n",
                    _qname);
            exit(1);
        }

        sprintf(line, "%Zu\t%Zu\t%s\t%Zu\t%Zu\t%s\t%Zu",
                fragment_id, _flagval, _rname, reported_pos,
                _mapq, _cigar, _tagval);
    }
    else
    {
        // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
        char rnext[] = "=";
        size_t pnext = 0;
        int tlen = 0;
        sprintf(line, "%s\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i\t%s\t%s",
                _qname, flag.get_raw(), _rname, reported_pos,
                _mapq, _cigar, rnext, pnext, tlen, "*", "*");
    }

    if (strlen(_tag_string) != 0)
    {
        strcat(line, "\t");
        strcat(line, _tag_string);
    }

    this->Init(line);
}


SamLine::SamLine(FILE * seqfile)
{
    char buf[4096 + 1];
    if (fgets(buf, 4096, seqfile) == NULL)
    {
        this->parse_flag = END_OF_FILE;
        return;
    }
    char * nl = strchr(buf, '\n');
    assert(nl != NULL);
    *nl = '\0';

    this->Init(buf);
}



SamLine::SamLine(char const* samline_string)
{
    this->Init(samline_string);
}



SamLine::SamLine(SamLine const& s) :
    parse_flag(s.parse_flag),
    fragment_id(s.fragment_id),
    flag(s.flag),
    rname(NULL),
    pos(s.pos),
    mapq(s.mapq),
    cigar(NULL),
    rnext(NULL),
    pnext(s.pnext),
    tlen(s.tlen),
    tags(s.tags),
    tag_string(NULL),
    flattened_pos(s.flattened_pos),
    cigar_compared(NULL),
    qname_string(NULL)
{
    if (s.rname != NULL)
    {
        this->rname = new char[strlen(s.rname) + 1];
        strcpy(this->rname, s.rname);
    }
    if (s.cigar != NULL)
    {
        this->cigar = new char[strlen(s.cigar) + 1];
        strcpy(this->cigar, s.cigar);
    }
    if (s.cigar_compared != NULL)
    {
        this->cigar_compared = new char[strlen(s.cigar_compared) + 1];
        strcpy(this->cigar_compared, s.cigar_compared);
    }
    if (s.qname_string != NULL)
    {
        this->qname_string = new char[strlen(s.qname_string) + 1];
        strcpy(this->qname_string, s.qname_string);
    }
    if (s.rnext != NULL)
    {
        this->rnext = new char[strlen(s.rnext) + 1];
        strcpy(this->rnext, s.rnext);
    }
    if (s.tag_string != NULL)
    {
        this->tag_string = new char[strlen(s.tag_string) + 1];
        strcpy(this->tag_string, s.tag_string);
    }
}


//Merge two or more SAMlines that are part of the same fragment
SamLine::SamLine(SamLine const* samlines[], size_t n_lines, char const* read_layout)
    : parse_flag(samlines[0]->parse_flag),
      fragment_id(samlines[0]->fragment_id),
      pos(samlines[0]->pos),
      mapq(samlines[0]->mapq),
      rnext(NULL),
      pnext(0),
      tlen(0),
      flattened_pos(samlines[0]->flattened_pos),
      cigar_compared(NULL),
      qname_string(NULL) 
// since we are merging SAM (not rSAM), we don't yet have a use for
// cigar_compared
{
    this->flag.set_raw(0);

    //assume at least one line
    assert(n_lines > 0);

    SamLine const* sl = samlines[0];

    this->tags = sl->tags;


    //assume valid samlines, and that they are sorted by alignment position.

    // b. construct merged CIGAR
    char merged_cigar[1024];

    char cigar_op[16];
    Cigar::CIGAR_VEC prev_cigar = Cigar::FromString(sl->cigar, 0);

    strcpy(merged_cigar, Cigar::ToString(prev_cigar).c_str());

    Cigar::CIGAR_VEC cur_cigar;

    // take an 'innocent until proven guilty strategy here
    // 
    this->flag.all_fragments_mapped = (! sl->flag.this_fragment_unmapped);

    for (size_t i = 1; i != n_lines; ++i)
    {
        this->tags.raw_score += samlines[i]->tags.raw_score;

        int64_t inter_seq_jump = 
            static_cast<int64_t>(samlines[i]->pos) 
            - (static_cast<int64_t>(samlines[i-1]->pos) 
               + static_cast<int64_t>(Cigar::Length(prev_cigar, true)));

        cur_cigar = Cigar::FromString(samlines[i]->cigar, 0);
        
        sprintf(cigar_op, "%ziT", inter_seq_jump);
        strcat(merged_cigar, cigar_op);
        strcat(merged_cigar, Cigar::ToString(cur_cigar).c_str());

        prev_cigar = cur_cigar;

        this->flag.all_fragments_mapped = 
            this->flag.all_fragments_mapped &&
            (! samlines[i]->flag.this_fragment_unmapped);
    }
    // c. merge tags appropriately (complicated!)  perhaps for now, we
    // simply can catenate tag strings together in fragment alignment
    // order. then, when printing, de-duplicate them...

    // for now, just take the tag string of the first record as the
    // string for the whole one
    
    //2 initialize all remaining fields and pointers

    this->rname = new char[strlen(sl->rname) + 1];
    strcpy(this->rname, sl->rname);

    this->cigar = new char[strlen(merged_cigar) + 1];
    strcpy(this->cigar, merged_cigar);

    this->tag_string = new char[strlen(sl->tag_string) + 1];
    strcpy(this->tag_string, sl->tag_string);

    this->flag.is_rsam_format = 1;
    this->flag.num_fragments_in_template = n_lines;

    size_t fi = sl->flag.first_fragment_in_template ? 0 : n_lines - 1;
    int sl_frag_to_temp_layout = (read_layout[fi] == 'f') ?
        LAYOUT_FORWARD : LAYOUT_REVERSE;

    int sl_frag_to_ref_layout = 
        sl->flag.this_fragment_on_neg_strand ?
        LAYOUT_REVERSE : LAYOUT_FORWARD;

    this->flag.template_layout = 
        (sl_frag_to_temp_layout == sl_frag_to_ref_layout) 
        ? LAYOUT_FORWARD
        : LAYOUT_REVERSE;

    for (size_t i = 0; i != n_lines; ++i)
    {
        this->flag.read_layout |=
            (read_layout[i] == 'f') ? (LAYOUT_FORWARD<<i) : (LAYOUT_REVERSE<<i);
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
    this->flag.alignment_not_primary = sl->flag.alignment_not_primary;

    // these fields are not applicable to rSAM format
    this->flag.multi_fragment_template = 0;
    this->flag.this_fragment_unmapped = 0;
    this->flag.next_fragment_unmapped = 0;
    this->flag.this_fragment_on_neg_strand = 0;
    this->flag.next_fragment_on_neg_strand = 0;
    this->flag.first_fragment_in_template = 0;
    this->flag.last_fragment_in_template = 0;

}


// must be null-terminated
void SamLine::Init(char const* samline_string)
{

    assert(strlen(samline_string) < 4096);
    assert(strchr(samline_string, '\n') == NULL); // must not contain a newline
    this->qname_string = NULL;

    if (samline_string[0] == '@')
    {
        //this is a header
        this->parse_flag = HEADER;

        this->tag_string = new char[strlen(samline_string) + 1];
        strcpy(this->tag_string, samline_string + 1);
        return;
    }

    else
    {
        bool was_error = false;
        int tag_start;

        if (SamLine::expect_rsam_format)
        {
            //FID.FLAG.RNAME.POS.MAPQ.CIGAR.TAGRAW[.TAG[.TAG[.TAG...]]]
            size_t tag_raw;
            size_t flag_raw;

            int num_fields =
                sscanf(samline_string, 
                       "%zu\t" "%zu\t" "%a[^\t]\t" "%zu\t" "%zu\t"
                       "%as\t" "%zu%n",
                       &this->fragment_id, &flag_raw, &this->rname, &this->pos, &this->mapq, 
                       &this->cigar, &tag_raw, &tag_start);

            this->flag.set_raw(flag_raw);
            this->tags.set_raw(tag_raw);

            // samline_string + tag_start either is a '\t' or NULL
            char const* tag_ptr = strlen(samline_string + tag_start) > 0
                ? samline_string + tag_start + 1
                : NULL;

            if (tag_ptr != NULL)
            {
                this->tag_string = new char[strlen(tag_ptr) + 1];
                strcpy(this->tag_string, tag_ptr);
            }
            else
            {
                this->tag_string = NULL;
            }

            if (num_fields != 7)
            {
                was_error = true;
            }
            else
            {
                this->rnext = NULL;
                this->pnext = 0;
            }
        }
        else
        {
            //QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
            char * qname;
            int seq_start, seq_end;
            size_t flag_raw;

            int num_fields =
                sscanf(samline_string, 
                       "%a[^\t]\t" "%zu\t" "%a[^\t]\t" "%zu\t" "%zu\t" 
                       "%as\t" "%as\t" "%zu\t" "%i\t"
                       "%n%*s%n\t%*s" 
                       "%n",
                       &qname, &flag_raw, &this->rname, &this->pos, &this->mapq, 
                       &this->cigar, &this->rnext, &this->pnext, &this->tlen,
                       &seq_start, &seq_end,
                       &tag_start);

            this->flag.set_raw(flag_raw);

            // initialize tags and tag_string
            init_tags_from_string(samline_string + tag_start, 
                                  SamLine::raw_score_tag,
                                  SamLine::worst_fragment_score,
                                  & this->tags,
                                  & this->tag_string);
        
            if (num_fields != 9)
            {
                was_error = true;
            }
            else
            {
                this->fragment_id = this->parse_fragment_id(qname);

                if (this->fragment_id == QNAME_FORMAT_ERROR)
                {
                    fprintf(stderr, "Error: samline_fragment_id: bad qname format in this line\n%s\n\n",
                            qname);
                    exit(1);
                }

                if (SamLine::retain_qname_string)
                {
                    this->qname_string = qname;
                }
                else
                {
                    this->qname_string = NULL;
                    delete [] qname;
                }
                if (strcmp(this->cigar, "*") == 0)
                {
                    //replace the '*' CIGAR with the more useful 'nS' op.
                    assert(this->flag.this_fragment_unmapped);

                    char cigar_tmp[32];
                    sprintf(cigar_tmp, "%iS", seq_end - seq_start);

                    delete [] this->cigar;
                    this->cigar = new char[strlen(cigar_tmp) + 1];
                    strcpy(this->cigar, cigar_tmp);
                }
            }
        }

        if (was_error)
        {
            char const* format = SamLine::expect_rsam_format ? "rSAM" : "SAM";
            fprintf(stderr, "Error: SamLine: Bad %s format:\n\n%s\n",
                    format, samline_string);
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

        this->parse_flag = DATA_LINE;

        this->flattened_pos = 0;
        this->cigar_compared = NULL;
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


// splice out zero-length 
void SamLine::SetCigarCompared()
{
    assert(this->cigar_compared == NULL);

    // convert 'T' to 'N' and consolidate.  T and N are the only ops
    // that consume ref but not seq.  They differ in whether template
    // is consumed.  By collapsing T to N, two CIGARs that differ in
    // the alignment of their unsequenced portion only, will be made
    // identical.
    size_t cigar_strlen = strlen(this->cigar);
    char * cigar_in = new char[cigar_strlen + 1];
    char * cigar_out = new char[cigar_strlen + 1];
    strcpy(cigar_in, this->cigar);

    char * c = cigar_in;
    while ((c = strchr(c, 'T')) != NULL)
    {
        *c++ = 'N';
    }
    Cigar::Consolidate(cigar_in, cigar_out);
    this->cigar_compared = new char[strlen(cigar_out) + 1];
    strcpy(this->cigar_compared, cigar_out);
    delete cigar_in;
    delete cigar_out;
}

SamLine::~SamLine()
{
    char ** fields[] = { & this->rname, & this->cigar, 
                         & this->rnext, & this->tag_string, 
                         & this->cigar_compared, & this->qname_string, 
                         NULL };

    for (char *** f = fields; *f != NULL; ++f)
    {
        if (**f != NULL)
        {
            delete [] **f;
            **f = NULL;
        }
    }
}


void SamLine::fprint(FILE * sam_fh) const
{
    this->print_aux(sam_fh, true);
}

void SamLine::sprint(char * buf) const
{
    this->print_aux(buf, false);
}

void SamLine::print_aux(void * print_buf, bool to_file) const
{

    if (this->parse_flag == HEADER)
    {
        FileUtils::switch_printf(to_file, & print_buf, "%s\n", this->tag_string);
    }

    else if (this->flag.is_rsam_format)
    {
        // is rSAM.
        size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();
        
        // FID.FLAG.RNAME.POS.MAPQ.CIGAR[.TAG[.TAG[.TAG...]]]
        FileUtils::switch_printf(to_file, & print_buf, "%zu\t%zu\t%s\t%zu\t%zu\t%s\t%zu",
                                 this->fragment_id, this->flag.get_raw(), this->rname, reported_pos,
                                 this->mapq, this->cigar, this->tags.get_raw());
        
        if (this->tag_string != NULL)
        {
            FileUtils::switch_printf(to_file, & print_buf, "\t%s", this->tag_string);
        }
        FileUtils::switch_printf(to_file, & print_buf, "\n");
    }
    else
    {
        // is traditional SAM.
        size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();
        size_t reported_pnext = this->flag.next_fragment_unmapped ? 0 : this->ones_based_pnext();

        // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
        FileUtils::switch_printf(to_file, & print_buf, 
                                 "%zu\t%zu\t%s\t%zu\t%zu\t%s\t%s\t%zu\t%i\t%s\t%s",
                                 this->fragment_id, flag.get_raw(), this->rname, reported_pos,
                                 this->mapq, this->cigar, this->rnext, reported_pnext, this->tlen,
                                 ".", ".");
        
        if (this->tag_string != NULL)
        {
            FileUtils::switch_printf(to_file, & print_buf, "\t%s", this->tag_string);
        }
        
        FileUtils::switch_printf(to_file, & print_buf, "\n");
    }
    return;
}




// print an rSAM record as a set of SAM records, to be stored in
// 'out_string'. 'seq_data' holds tab-separated fields of: qname seq
// qual [seq qual...]  and must contain as many seq/qual pairs as
// there are fragments on the template Will print appropriately
// reverse complemented SEQ / QUAL fields, depending on the flags.
void SamLine::print_rsam_as_sam(char const* seq_data, char * out_string) const
{
    int64_t cur_read_start = this->pos;

    char rnext[2] = "=";
    
    Cigar::CIGAR_VEC cigar_v = Cigar::FromString(this->cigar, 0);

    Cigar::CIGAR_ITER read_start = cigar_v.begin();
    Cigar::CIGAR_ITER read_end = cigar_v.begin();
    Cigar::CIGAR_ITER unseq_end = cigar_v.begin();
    
    int64_t tlen = Cigar::Length(cigar_v, false);

    char cigar_substring[1024];

    SamFlag outflag = this->flag;

    outflag.template_layout = 0;
    outflag.is_rsam_format = 0;
    outflag.num_fragments_in_template = 0;
    outflag.read_layout = 0;

    size_t num_fragments = this->flag.num_fragments_in_template;

    outflag.multi_fragment_template = (num_fragments > 1) ? 1 : 0;

    // find all seqs and quals in 'seq_data'
    char const** seqs = new char const*[num_fragments];
    char const** quals = new char const*[num_fragments];
    char const* qname = seq_data;

    // point seqs and quals into the seq_data array
    char const* ptr = seq_data;
    for (size_t f = 0; f != num_fragments * 2; ++f)
    {
        ptr = strchr(ptr, '\t');
        if (ptr == NULL)
        {
            fprintf(stderr, "Error: print_rsam_as_sam: there are %zu fragments but\n"
                    "seq_data should have %zu fields (qname, seq, qual [, seq, qual ...]\n"
                    "bad seq_data line is\n%s\n",
                    num_fragments, 1 + (num_fragments * 2), seq_data);
            exit(1);
        }
        else
        {
            ++ptr;
        }
        (f % 2 == 0) ? seqs[f/2] = ptr : quals[f/2] = ptr;
    }

    size_t align_index = 0; // index in left-to-right alignment order
    size_t layout_index;
    int this_fragment_layout;
    int next_fragment_layout;

    size_t reported_pos;
    size_t reported_pnext;

    while (read_start != cigar_v.end())
    {
        layout_index = (this->flag.template_layout == LAYOUT_FORWARD) 
            ? align_index 
            : num_fragments - align_index - 1;

        this_fragment_layout = 1 & (this->flag.read_layout>>layout_index);
        next_fragment_layout = 1 & (this->flag.read_layout>>((layout_index + 1) % num_fragments));

        //construct the next [read_start, read_end) sequenced portion range
        assert((*read_start).op.seq); // want to be within the sequenced portion
        while (read_end != cigar_v.end() && (*read_end).op.temp == (*read_end).op.seq)
        {
            ++read_end;
        }
        unseq_end = read_end;
        while (unseq_end != cigar_v.end() && ! (*unseq_end).op.seq)
        {
            ++unseq_end;
        }
        // now [read_start, read_end) is a read alignment range
        // and [read_end, unseq_end) is the unsequenced portion alignment range
        size_t start_to_start_dist = Cigar::Length(read_start, unseq_end, true);
        if (this->flag.all_fragments_mapped)
        {
            Cigar::Trim(read_start, read_end, false, & read_start, & read_end);
            Cigar::ToString(read_start, read_end, cigar_substring);
            reported_pos = cur_read_start + 1;
            reported_pnext = (read_end == cigar_v.end())
                ? this->pos + 1
                : this->pos + start_to_start_dist + 1;
        }
        else
        {
            strcpy(cigar_substring, "*");
            reported_pos = 0;
            reported_pnext = 0;
        }

        outflag.this_fragment_unmapped = 1 - this->flag.all_fragments_mapped;
        outflag.next_fragment_unmapped = 1 - this->flag.all_fragments_mapped;

        //fragment_layout is fragment-to-template.
        //template_layout is template-to-reference.
        //we want to calculate fragment-to-reference.
        outflag.this_fragment_on_neg_strand =
            (this_fragment_layout == this->flag.template_layout) 
            ? LAYOUT_FORWARD
            : LAYOUT_REVERSE;

        outflag.next_fragment_on_neg_strand =
            (next_fragment_layout == this->flag.template_layout) 
            ? LAYOUT_FORWARD
            : LAYOUT_REVERSE;

        outflag.first_fragment_in_template = (layout_index == 0);
        outflag.last_fragment_in_template = (layout_index == num_fragments - 1);
        
        // definition of tlen in SAM 1.4 Spec
        int64_t used_tlen;
        if ((! this->flag.all_fragments_mapped)
            || num_fragments == 1 
            || (layout_index != 0 && layout_index != (num_fragments - 1)))
        {
            used_tlen = 0;
        }
        else
        {
            used_tlen = (align_index == 0) ? tlen : - tlen;
        }

        char const* qname_end = strchr(qname, '\t');
        std::copy(qname, qname_end, out_string);
        out_string += std::distance(qname, qname_end);
        out_string += 
            sprintf(out_string, 
                    "\t%zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%zi",
                    outflag.get_raw(), this->rname, reported_pos,
                    this->mapq, cigar_substring, rnext, reported_pnext, used_tlen);
        
        size_t seq_len = std::distance(seqs[layout_index], quals[layout_index]) - 1;
        if (outflag.this_fragment_on_neg_strand)
        {
            *out_string++ = '\t';
            reverse_comp(seqs[layout_index], seqs[layout_index] + seq_len, out_string);
            out_string += seq_len;
            *out_string++ = '\t';
            std::reverse_copy(quals[layout_index], quals[layout_index] + seq_len, out_string);
            out_string += seq_len;
        }
        else
        {
            *out_string++ = '\t';
            std::copy(seqs[layout_index], seqs[layout_index] + seq_len, out_string);
            out_string += seq_len;
            *out_string++ = '\t';
            std::copy(quals[layout_index], quals[layout_index] + seq_len, out_string);
            out_string += seq_len;
        }
        
        //now about the tags.
        // this should be optional, 
        if (this->tags.raw_score_present)
        {
            out_string += sprintf(out_string, "\t%c%c:i:%i", 
                                  this->tags.raw_score_tag[0], 
                                  this->tags.raw_score_tag[1],
                                  this->tags.raw_score);
        }
        if (this->tags.alignment_space_present)
        {
            out_string += sprintf(out_string, "\t%s:A:%c", 
                                  AlignSpaceTag, this->tags.alignment_space);
        }
        if (this->tags.stratum_rank_present)
        {
            out_string += sprintf(out_string, "\t%s:i:%i", 
                                  StratumRankTag, this->tags.stratum_rank);
        }
        if (this->tags.stratum_size_present)
        {
            out_string += sprintf(out_string, "\t%s:i:%i", 
                                  StratumSizeTag, this->tags.stratum_size);
        }
        if (this->tag_string != NULL)
        {
            out_string += sprintf(out_string, "\t%s", this->tag_string);
        }

        out_string += sprintf(out_string, "\n");

        cur_read_start += start_to_start_dist;
        read_start = unseq_end;
        read_end = unseq_end;
        ++align_index;
    }
    delete [] seqs;
    delete [] quals;
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

        delete [] this->tag_string;
        this->tag_string = new_extra_tag;
    }
}


size_t SamLine::template_length() const 
{ 
    return Cigar::Length(Cigar::FromString(this->cigar, 0), false); 
}


//true if samlines represent sequence mate pairs, regardless of
//whether they are mapped as such
bool AreSequencedMatePairs(SamLine const& a, SamLine const& b)
{
    assert(! (a.flag.is_rsam_format || b.flag.is_rsam_format));

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


// reads and allocates a buffer loaded with the sam header
char * ReadAllocSAMHeader(FILE * sam_fh)
{
    std::vector<char> header;
    header.reserve(1000);

    char p,n;
    p = '\0';

    while (! feof(sam_fh))
    {
        n = fgetc(sam_fh);
        if ((p == '\0' || p == '\n') && n != '@')
        {
            ungetc(n, sam_fh);
            clearerr(sam_fh);
            break;
        }
        p = n;
        header.push_back(p);
    }
    char * header_string = new char[header.size() + 1];
    std::copy(header.begin(), header.end(), header_string);
    header_string[header.size()] = '\0';

    return header_string;
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


char const* SamLine::cigar_for_comparison() const
{
    return this->cigar_compared == NULL
        ? this->cigar
        : this->cigar_compared;
}

void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh)
{
    SetToFirstDataLine(input_sam_fh);

    size_t header_length = ftell(*input_sam_fh);
    rewind(*input_sam_fh);

    char * header_buf = new char[header_length];
    fread(header_buf, 1, header_length, *input_sam_fh);
    fwrite(header_buf, 1, header_length, output_fh);
    delete [] header_buf;
    fflush(output_fh);
}



