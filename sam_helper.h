#ifndef _SAM_HELPER_H
#define _SAM_HELPER_H

#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <cassert>
#include <map>
#include <functional>

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

//flag masks.  When set, these flags proclaim their name
/* namespace SamFlags */
/* { */
/*     extern int const MULTI_FRAGMENT_TEMPLATE; */
/*     extern int const ALL_FRAGMENTS_MAPPED; */
/*     extern int const THIS_FRAGMENT_UNMAPPED; */
/*     extern int const NEXT_FRAGMENT_UNMAPPED; */
/*     extern int const THIS_FRAGMENT_ON_NEG_STRAND; */
/*     extern int const NEXT_FRAGMENT_ON_NEG_STRAND; */
/*     extern int const FIRST_FRAGMENT_IN_TEMPLATE; */
/*     extern int const LAST_FRAGMENT_IN_TEMPLATE; */
/*     extern int const ALIGNMENT_NOT_PRIMARY; */
/*     extern int const FAILED_QUALITY_CHECK; */
/*     extern int const PCR_OR_OPTICAL_DUPLICATE; */
    
/*     extern int const TEMPLATE_ON_NEG_STRAND; */
/* }; */

enum SAM_PARSE
    {
        END_OF_FILE,
        PARSE_ERROR,
        HEADER,
        DATA_LINE
    };

enum SAM_QNAME_FORMAT
    {
        SAM_NUMERIC,
        SAM_ILLUMINA,
        SAM_CASAVA18,
        SAM_NON_INTERPRETED
    };

struct less_char_ptr
{
    bool operator()(char const* a, char const* b) const
    {
        return strcmp(a, b) < 0;
    }
};


struct eqstr
{
    bool operator()(const char* s1, const char* s2) const;
};


#ifdef __GXX_EXPERIMENTAL_CXX0X__
typedef std::hash<std::string> STRING_HASH;
struct to_integer : public STRING_HASH
{
    size_t operator()(char const* k) const;
};
typedef std::unordered_map<char const*, char, to_integer, eqstr> CONTIG_SPACE;
typedef std::unordered_map<char const*, size_t, to_integer, eqstr> CONTIG_OFFSETS;
#else
typedef std::tr1::hash<std::string> STRING_HASH;
struct to_integer : public STRING_HASH
{
    size_t operator()(char const* k) const;
};
typedef std::tr1::unordered_map<char const*, char, to_integer, eqstr> CONTIG_SPACE;
typedef std::tr1::unordered_map<char const*, size_t, to_integer, eqstr> CONTIG_OFFSETS;
#endif

typedef CONTIG_OFFSETS::const_iterator OFFSETS_ITER;



#define AlignSpaceTag "XP"
#define AlignSpaceType 'A'
#define AlignSpaceMissing '-'

// The Not Applicable Read Layout used when parsing 
#define ReadLayoutNA "-"

template<typename V>
V parse_sam_tag(char const* tag_string, 
                char const* tag_name, 
                char const* tag_format, 
                V default_if_missing, bool * present)
{
    V tag_value;
    char const* tag_substring = tag_string == NULL ? NULL : strstr(tag_string, tag_name);

    char tag_plus_format[32];
    strcpy(tag_plus_format, tag_name);
    strcat(tag_plus_format, tag_format);

    if (tag_substring != NULL 
        && sscanf(tag_substring, tag_plus_format, &tag_value) == 1)
    {
        *present = true;
        return tag_value;
    }
    else
    {
        *present = false;
        return default_if_missing;
    }
}


/*
  Storage policy:
  line, extra, and extra_tag own the memory they point to.
  all other pointers point to offsets within them.
  qname, rname, rnext, seq, and qual always point into line.
  cigar may point into line or extra
  tag_string may point into line or extra_tag

  SAM Format (1.4-r983) (tab-separated, including optional tags, represented here by '.')
  QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN.SEQ.QUAL[.TAG[.TAG[.TAG...]]]

  Or, with zero-length SEQ and QUAL:
  QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN..[.TAG[.TAG[.TAG...]]]

  Or, with no tags and zero length SEQ and QUAL:
  QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.RNEXT.PNEXT.TLEN..

  rSAM Format (recombinant, template-based SAM)
  QNAME.FLAG.RNAME.POS.MAPQ.CIGAR.SEQ.QUAL[.TAG[.TAG[.TAG...]]]
  
 */

union SamFlag
{
    struct
    {
        unsigned int multi_fragment_template : 1;
        unsigned int all_fragments_mapped : 1;
        unsigned int this_fragment_unmapped : 1;
        unsigned int next_fragment_unmapped : 1;
        unsigned int this_fragment_on_neg_strand : 1;
        unsigned int next_fragment_on_neg_strand : 1;
        unsigned int first_fragment_in_template : 1;
        unsigned int last_fragment_in_template : 1;
        unsigned int alignment_not_primary : 1;
        unsigned int failed_quality_check : 1;
        unsigned int pcr_or_optical_duplicate : 1;
        unsigned int template_layout : 1;
        unsigned int : 4; // padding to 16 bits

        unsigned int is_rsam_format : 8;
        unsigned int num_fragments_in_template : 8;
        unsigned int read_layout : 32;
    };
    size_t raw;

};

#define LAYOUT_FORWARD 0
#define LAYOUT_REVERSE 1




class SamLine
{
public:

    size_t bytes_in_line;
    char * line;
    SAM_PARSE parse_flag;
    char * qname;
    size_t fragment_id;
    SamFlag flag;
    char * rname;
    size_t pos;
    size_t mapq;
    char * cigar;
    char * rnext;
    size_t pnext;
    int tlen;
    char * seq;
    char * qual;
    char * tag_string;
    char * extra;
    char * extra_tag;
    size_t flattened_pos; //position along a virtual concatenated meta-contig.
    char alignment_space;
    

    SamLine(SAM_PARSE _parse_flag,
            char const* _qname, size_t _flag, 
            char const* _rname, size_t _pos,
            size_t _mapq, char const* cigar,
            char const* _seq,
            char const* _qual,
            char const* _tag_string);

    SamLine(FILE * sam_fh);
    SamLine(char const* samline_string);

    SamLine(SamLine const& s);

    SamLine(SamLine const* samlines[], size_t n_lines, char const* read_layout);

    void Init(char const* samline_string);
    void SetFlattenedPosition(CONTIG_OFFSETS const& contig_offsets,
                              CONTIG_OFFSETS::const_iterator * contig_iter);

    static void SetGlobalFlags(SAM_QNAME_FORMAT _qname_format,
                               char const* _expected_layout);

    static SAM_QNAME_FORMAT sam_qname_format;
    static char expected_read_layout[256];
    static bool expect_rsam_format;

    static size_t (* parse_fragment_id)(char const* qname);

    SamLine();

    ~SamLine();


    //bool operator<(SamLine const& b) const;
    //bool operator==(SamLine const& b) const;

    SAM_PARSE load(FILE * seqfile);

    void print(FILE * samfile, bool print_rsam) const;
    void print_sam(FILE * samfile) const;
    void print_rsam(FILE * samfile) const;

    void print_fastq(FILE * fastq_fh) const;

    int alignment_score(char const* tag, int default_if_missing, bool * has_score) const;
    //char alignment_space(char default_if_missing, bool * has_space) const;

    //queries the samline for the desired tag code, loading type and value_string if found.
    //returns whether the tag was found
    bool has_tag(char const* tag_code, char & type, char * value_string) const;

    //adds the desired tag. 
    void add_tag(char const* tag_code, char type, char const* value_string);

    size_t aligned_read_length();
    size_t raw_read_length() const { return strlen(this->seq); }
    size_t ones_based_pos() const { return this->pos + 1; }
    size_t zero_based_pos() const { return this->pos; }

    /* bool multi_fragment_template() const { return (this->flag & SamFlags::MULTI_FRAGMENT_TEMPLATE) != 0; } */
    /* bool all_fragments_mapped() const { return (this->flag & SamFlags::ALL_FRAGMENTS_MAPPED) != 0; } */
    /* bool this_fragment_unmapped() const { return ((this->flag & SamFlags::THIS_FRAGMENT_UNMAPPED)) != 0; } */
    /* bool next_fragment_unmapped() const { return (this->flag & SamFlags::NEXT_FRAGMENT_UNMAPPED) != 0; } */

    /* //strandedness flags are nonzero for the negative strand. */
    /* bool template_on_pos_strand() const { return (this->flag & SamFlags::TEMPLATE_ON_NEG_STRAND) == 0; } */
    /* bool this_fragment_on_pos_strand() const { return (this->flag & SamFlags::THIS_FRAGMENT_ON_NEG_STRAND) == 0; } */
    /* bool next_fragment_on_pos_strand() const { return (this->flag & SamFlags::NEXT_FRAGMENT_ON_NEG_STRAND) == 0; } */

    /* bool first_fragment_in_template() const { return (this->flag & SamFlags::FIRST_FRAGMENT_IN_TEMPLATE) != 0; } */
    /* bool last_fragment_in_template() const { return (this->flag & SamFlags::LAST_FRAGMENT_IN_TEMPLATE) != 0; } */
    /* bool alignment_not_primary() const { return (this->flag & SamFlags::ALIGNMENT_NOT_PRIMARY) != 0; } */
    /* bool failed_quality_check() const { return (this->flag & SamFlags::FAILED_QUALITY_CHECK) != 0; } */
    /* bool pcr_or_optical_duplicate() const { return (this->flag & SamFlags::PCR_OR_OPTICAL_DUPLICATE) != 0; } */

    char const* next_fragment_ref_name() const;

};


void SetToFirstDataLine(FILE ** sam_fh);




enum Strand
    {
        POS,
        NEG,
        POS_NEG
    };


//test whether a and b are indeed mate pairs
bool AreSequencedMatePairs(SamLine const& a, SamLine const& b);

bool AreMappedMatePairs(SamLine const& a, SamLine const& b);

bool AreUnmappedMatePairs(SamLine const& a, SamLine const& b);


char const* strand_to_char(Strand strand);

int SAM_cmp_qname_flag(SamLine const& a, SamLine const& b);

int SAM_cmp_qname_flag_aux(char const* qname1, int flag1,
                           char const* qname2, int flag2);


void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh);

#endif // _SAM_HELPER_H

