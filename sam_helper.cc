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


void reverse_comp_inplace(char const* begin, char const* end, char * rcomp)
{

    while (begin != end)
    {
        --end;
        *rcomp = base_to_complement[static_cast<int>(*end)];
        ++rcomp;
    }
}


void encode_tags(char const* tag_string, 
                 char const* raw_score_tag,
                 size_t worst_raw_score, 
                 SamTag * tags,
                 char ** extra_tags)
{
    char tmp_buf[1024];
    char * tmp_buf_ptr = tmp_buf;
    tmp_buf_ptr[0] = '\0';

    char const* tag = tag_string;
    char tag_name[3];
    char tag_type;
    char tag_value[1024];
    int advance;

    while (sscanf(tag, " %[^:]:%c:%s%n", tag_name, &tag_type, tag_value, &advance) == 3)
    {
        if (strcmp(tag_name, raw_score_tag) == 0)
        {
            (*tags).raw_score = atoi(tag_value);
        }
        else if (strcmp(tag_name, AlignSpaceTag) == 0)
        {
            (*tags).alignment_space = tag_value[0];
        }
        else if (strcmp(tag_name, StratumRankTag) == 0)
        {
            (*tags).stratum_rank = atoi(tag_value);
        }
        else if (strcmp(tag_name, StratumSizeTag) == 0)
        {
            (*tags).stratum_size = atoi(tag_value);
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

size_t (* SamLine::parse_fragment_id)(char const* qname) = &parse_fragment_id_zero;



void SamLine::SetGlobalFlags(SAM_QNAME_FORMAT _qname_format,
                             char const* _expected_layout,
                             char const* _raw_score_tag,
                             size_t _worst_fragment_score)
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
    flag.raw = _flagval;

    size_t reported_pos = flag.this_fragment_unmapped ? 0 : _pos + 1;

    if (SamLine::expect_rsam_format)
    {
        // FID.FLAG.RNAME.POS.MAPQ.CIGAR.TAGRAW[.TAG[.TAG[.TAG...]]]
        size_t fragment_id = this->parse_fragment_id(_qname);

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
                _qname, flag.raw, _rname, reported_pos,
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
    flattened_pos(s.flattened_pos)
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
      flattened_pos(samlines[0]->flattened_pos)
{
    this->flag.raw = 0;

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

    for (size_t i = 1; i != n_lines; ++i)
    {
        assert(samlines[i]->tags.alignment_space == this->tags.alignment_space);
        assert(samlines[i]->tags.stratum_rank == this->tags.stratum_rank);
        assert(samlines[i]->tags.stratum_size == this->tags.stratum_size);

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

    this->flag.raw = 0;

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
    this->flag.all_fragments_mapped = sl->flag.all_fragments_mapped;
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
            int num_fields =
                sscanf(samline_string, 
                       "%zu\t" "%zu\t" "%as\t" "%zu\t" "%zu\t"
                       "%as\t" "%zu%n",
                       &this->fragment_id, &this->flag.raw, &this->rname, &this->pos, &this->mapq, 
                       &this->cigar, &this->tags.raw, &tag_start);

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

            int num_fields =
                sscanf(samline_string, 
                       "%as\t" "%zu\t" "%as\t" "%zu\t" "%zu\t" 
                       "%as\t" "%as\t" "%zu\t" "%i\t"
                       "%n%*s%n\t%*s" 
                       "%n",
                       &qname, &this->flag.raw, &this->rname, &this->pos, &this->mapq, 
                       &this->cigar, &this->rnext, &this->pnext, &this->tlen,
                       &seq_start, &seq_end,
                       &tag_start);

            this->tags.raw_score = SamLine::worst_fragment_score;
            this->tags.alignment_space = AlignSpaceMissing;
            this->tags.stratum_rank = 0;
            this->tags.stratum_size = 0;

            if (num_fields != 9)
            {
                was_error = true;
            }
            else
            {
                this->fragment_id = this->parse_fragment_id(qname);

                delete [] qname;
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

        // initialize tag_string
        encode_tags(samline_string + tag_start, 
                    SamLine::raw_score_tag,
                    SamLine::worst_fragment_score,
                    & this->tags,
                    & this->tag_string);
        
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
    char ** fields[] = { & this->rname, & this->cigar, & this->rnext, & this->tag_string, NULL };

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
        FileUtils::switch_printf(to_file, print_buf, "%s\n", this->tag_string);
    }

    else if (this->flag.is_rsam_format)
    {
        size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();
        
        FileUtils::switch_printf(to_file, print_buf, "%zu\t%Zu\t%s\t%Zu\t%Zu\t%s\t%Zu",
                                 this->fragment_id, this->flag.raw, this->rname, reported_pos,
                                 this->mapq, this->cigar, this->tags.raw);
        
        if (this->tag_string != NULL)
        {
            FileUtils::switch_printf(to_file, print_buf, "\t%s", this->tag_string);
        }
        FileUtils::switch_printf(to_file, print_buf, "\n");
    }
    else
    {
        // is traditional SAM.
        size_t reported_pos = this->flag.this_fragment_unmapped ? 0 : this->ones_based_pos();
        size_t reported_pnext = this->flag.next_fragment_unmapped ? 0 : this->ones_based_pnext();

        // QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
        FileUtils::switch_printf(to_file, print_buf, 
                                 "%Zu\t%Zu\t%s\t%Zu\t%Zu\t%s\t%s\t%Zu\t%i\t%s\t%s",
                                 this->fragment_id, flag.raw, this->rname, reported_pos,
                                 this->mapq, this->cigar, this->rnext, reported_pnext, this->tlen,
                                 ".", ".");
        
        if (this->tag_string != NULL)
        {
            FileUtils::switch_printf(to_file, print_buf, "\t%s", this->tag_string);
        }
        
        FileUtils::switch_printf(to_file, print_buf, "\n");
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
    int64_t pnext;
    char rnext[2] = "=";
    
    Cigar::CIGAR_VEC cigar_v = Cigar::FromString(this->cigar, 0);

    Cigar::CIGAR_ITER start = cigar_v.begin();
    Cigar::CIGAR_ITER end = cigar_v.begin();

    int64_t tlen = Cigar::Length(cigar_v, false);

    char cigar_substring[1024];

    SamFlag flag = this->flag;
    flag.raw &= 255;

    size_t num_fragments = this->flag.num_fragments_in_template;

    flag.multi_fragment_template = (num_fragments > 1) ? 1 : 0;

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

    while (start != cigar_v.end())
    {
        layout_index = (this->flag.template_layout == LAYOUT_FORWARD) 
            ? align_index 
            : num_fragments - align_index - 1;

        this_fragment_layout = 1 & (this->flag.read_layout>>layout_index);
        next_fragment_layout = 1 & (this->flag.read_layout>>((layout_index + 1) % num_fragments));

        //construct the next [start, end) range for [^TU] CIGAR ops
        while (end != cigar_v.end() 
               && (*end).op.code != Cigar::T
               && (*end).op.code != Cigar::U)
        {
            ++end;
        }
        //now 'end' points to the real end or to a 'T' operator
        size_t read_length = Cigar::Length(start, end, false);

        if (this->flag.all_fragments_mapped)
        {
            Cigar::ToString(start, end, cigar_substring);
        }
        else
        {
            strcpy(cigar_substring, "*");
        }

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
        }

        flag.this_fragment_unmapped = 1 - this->flag.all_fragments_mapped;
        flag.next_fragment_unmapped = 1 - this->flag.all_fragments_mapped;

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

        flag.first_fragment_in_template = (layout_index == 0);
        flag.last_fragment_in_template = (layout_index == num_fragments - 1);
        
        // definition of tlen in SAM 1.4 Spec
        int64_t used_tlen;
        if (num_fragments == 1 
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
                    flag.raw, this->rname, cur_read_start + 1,
                    this->mapq, cigar_substring, rnext, pnext + 1, used_tlen);
        
        size_t seq_len = std::distance(seqs[layout_index], quals[layout_index]) - 1;
        if (flag.this_fragment_on_neg_strand)
        {
            *out_string++ = '\t';
            reverse_comp_inplace(seqs[layout_index], seqs[layout_index] + seq_len, out_string);
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
        out_string += 
            sprintf(out_string, "%s:i:%i\t%s:A:%c\t%s:i:%i\t%s:i:%i",
                    SamLine::raw_score_tag, this->tags.raw_score, 
                    AlignSpaceTag, this->tags.alignment_space,
                    StratumRankTag, this->tags.stratum_rank, 
                    StratumSizeTag, this->tags.stratum_size);

        if (this->tag_string != NULL)
        {
            sprintf(out_string, "\t%s", this->tag_string);
        }

        strcat(out_string, "\n");
        out_string += strlen(out_string);

        // update cur_read_start. [start, end) is a [^TU] range
        cur_read_start += Cigar::Length(start, end, true);
        start = end;
        while (end != cigar_v.end() 
               && ((*end).op.code == Cigar::T || (*end).op.code == Cigar::U))
        {
            ++end;
        }
        // update cur_read_start. [start, end) is a [TU] range
        cur_read_start += Cigar::Length(start, end, true);
        start = end;

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
    delete [] header_buf;
    fflush(output_fh);
}



