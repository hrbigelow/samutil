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

/*
  Policy: Assume either a SAM or rSAM format when parsing, according
  to the 'expect_rsam_format' flag.

  In either case, do not store qname, seq, or qual fields.
  Also, store 'unmapped' as 50S or 76S CIGAR (whatever the readlength is)
 */

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



// needed for buffer allocation estimate
size_t const MAX_RSAM_LINE_LENGTH = 20 + 20 + 20 + 12 + 3 + 100 + 20 + 12 + 300;

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

  rSAM Format (recombinant, template-based SAM)
  FID.FLAG.RNAME.POS.MAPQ.CIGAR[.TAG[.TAG[.TAG...]]]
  
 */

struct SamFlag
{
    unsigned int multi_fragment_template : 1;     // SAM
    unsigned int all_fragments_mapped : 1;        // SAM and rSAM
    unsigned int this_fragment_unmapped : 1;      // SAM
    unsigned int next_fragment_unmapped : 1;      // SAM
    unsigned int this_fragment_on_neg_strand : 1; // SAM
    unsigned int next_fragment_on_neg_strand : 1; // SAM
    unsigned int first_fragment_in_template : 1;  // SAM
    unsigned int last_fragment_in_template : 1;   // SAM
    unsigned int alignment_not_primary : 1;       // SAM and rSAM
    unsigned int failed_quality_check : 1;        // SAM and rSAM
    unsigned int pcr_or_optical_duplicate : 1;    // SAM and rSAM
    unsigned int template_layout : 1;             //         rSAM
    unsigned int : 4; // padding to 16 bits

    unsigned int is_rsam_format : 8;              //         rSAM
    unsigned int num_fragments_in_template : 8;   //         rSAM
    unsigned int read_layout : 32;                //         rSAM

    size_t get_raw() const;
    void set_raw(size_t raw);

};

struct SamTag
{
    uint16_t raw_score;
    uint16_t stratum_rank;
    uint16_t stratum_size;

    char alignment_space;

    bool raw_score_present;
    bool alignment_space_present;
    bool stratum_rank_present;
    bool stratum_size_present;

    char raw_score_tag[2];

    size_t get_raw() const;
    void set_raw(size_t raw);

};



#define AlignSpaceTag "XP"
#define StratumRankTag "XY"
#define StratumSizeTag "XZ"

#define AlignSpaceType 'A'
#define AlignSpaceMissing '-'
#define GuideAlignmentTag "XG"


#define LAYOUT_FORWARD 0
#define LAYOUT_REVERSE 1




class SamLine
{

    void print_aux(void * print_buf, bool to_file) const;

public:

    SAM_PARSE parse_flag;
    size_t fragment_id;
    SamFlag flag;
    char * rname;
    size_t pos;
    size_t mapq;
    char * cigar;
    char * rnext;
    size_t pnext;
    int tlen;
    SamTag tags;
    char * tag_string;
    size_t flattened_pos; //position along a virtual concatenated meta-contig.
    char * cigar_compared; // if not NULL, used for positional comparison

    SamLine(SAM_PARSE _parse_flag,
            char const* _qname, size_t _flag, 
            char const* _rname, size_t _pos,
            size_t _mapq, char const* cigar,
            size_t _tagval,
            char const* _tag_string);

    SamLine(FILE * sam_fh);
    SamLine(char const* samline_string);

    SamLine(SamLine const& s);

    SamLine(SamLine const* samlines[], size_t n_lines, char const* read_layout);

    void Init(char const* samline_string);
    void SetFlattenedPosition(CONTIG_OFFSETS const& contig_offsets,
                              CONTIG_OFFSETS::const_iterator * contig_iter);

    void SetCigarCompared();

    static void SetGlobalFlags(SAM_QNAME_FORMAT _qname_format,
                               char const* _expected_layout,
                               char const* _raw_score_tag,
                               size_t _worst_fragment_score);

    static SAM_QNAME_FORMAT sam_qname_format;
    static char expected_read_layout[256];
    static bool expect_rsam_format;
    static bool brief_records; // prints records with numeric qnames, and '*' for SEQ and QUAL
    static char raw_score_tag[3];
    static size_t worst_fragment_score;
    static bool initialized;

    static size_t (* parse_fragment_id)(char const* qname);

    SamLine();

    ~SamLine();


    //bool operator<(SamLine const& b) const;
    //bool operator==(SamLine const& b) const;

    SAM_PARSE load(FILE * seqfile);

    /* void print(FILE * samfile, bool print_rsam) const; */

    // print abbreviated SAM format to file
    void fprint(FILE * sam_fh) const;

    // print abbreviated SAM format to string
    void sprint(char * sam_string) const;

    void print_rsam_as_sam(char const* seq_data, char * out_string) const;

    void print_fastq(FILE * fastq_fh) const;

    int alignment_score(char const* tag, int default_if_missing, bool * has_score) const;
    //char alignment_space(char default_if_missing, bool * has_space) const;

    //queries the samline for the desired tag code, loading type and value_string if found.
    //returns whether the tag was found
    bool has_tag(char const* tag_code, char & type, char * value_string) const;

    //adds the desired tag. 
    void add_tag(char const* tag_code, char type, char const* value_string);

    size_t template_length() const;
    size_t ones_based_pos() const { return this->pos + 1; }
    size_t zero_based_pos() const { return this->pos; }

    size_t ones_based_pnext() const { return this->pnext + 1; }
    size_t zero_based_pnext() const { return this->pnext; }

    char const* next_fragment_ref_name() const;

    char const* cigar_for_comparison() const;
};


struct SamFilter
{
    size_t min_mapping_quality;
    size_t max_stratum_rank;
    size_t max_stratum_size;
    char const* alignment_space_filter;

    SamFlag comb;
    SamFlag values;

    SamFilter(char const* _tf,
              size_t _mmq,
              size_t _msr,
              size_t _mss,
              char const* _asf);

    bool pass(SamLine const* samline) const;
};




void SetToFirstDataLine(FILE ** sam_fh);

char * ReadAllocSAMHeader(FILE * sam_fh);



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

