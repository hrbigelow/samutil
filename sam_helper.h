#ifndef _SAM_HELPER_H
#define _SAM_HELPER_H

#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <cassert>
#include <map>


//#include <unordered_map>

//flag masks.  When set, these flags proclaim their name
namespace SamFlags
{
    extern int const PAIRED_IN_SEQUENCING;
    extern int const MAPPED_IN_PROPER_PAIR;
    extern int const QUERY_UNMAPPED;
    extern int const MATE_UNMAPPED;
    extern int const QUERY_ON_NEG_STRAND;
    extern int const MATE_ON_NEG_STRAND;
    extern int const FIRST_READ_IN_PAIR;
    extern int const SECOND_READ_IN_PAIR;
    extern int const ALIGNMENT_NOT_PRIMARY;
    extern int const FAILED_QUALITY_CHECK;
    extern int const PCR_OR_OPTICAL_DUPLICATE;
};

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



typedef std::map<char const*, size_t, less_char_ptr> CONTIG_OFFSETS;
typedef CONTIG_OFFSETS::const_iterator OFFSETS_ITER;


#define AlignSpaceTag "XP"
#define AlignSpaceType 'A'

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


//typedef std::unordered_map<char const*, size_t, std::hash<char const*>, eqstr> CONTIG_OFFSETS;

/*
  Storage policy:
  line, extra, and extra_tag own the memory they point to.
  all other pointers point to offsets within them.
  qname, rname, mrnm, seq, and qual always point into line.
  cigar may point into line or extra
  tag_string may point into line or extra_tag
 */

class SamLine
{
public:

    size_t bytes_in_line;
    char * line;
    SAM_PARSE parse_flag;
    char * qname;
    size_t fragment_id;
    int flag;
    char * rname;
    size_t pos;
    size_t mapq;
    char * cigar;
    char * mrnm;
    size_t mpos;
    int isize;
    char * seq;
    char * qual;
    char * tag_string;
    char * extra;
    char * extra_tag;
    size_t flattened_pos; //position along a virtual concatenated meta-contig.

    SamLine(SAM_PARSE _parse_flag,
            char const* _qname, int _flag, 
            char const* _rname, size_t _pos,
            size_t _mapq, char const* cigar,
            char const* _mrnm, size_t _mpos,
            int _isize, char const* _seq,
            char const* _qual,
            char const* _tag_string);

    SamLine(FILE * sam_fh, bool allow_absent_seq_qual);
    SamLine(char const* samline_string, bool allow_absent_seq_qual);

    SamLine(SamLine const& s);

    void Init(char const* samline_string, bool allow_absent_seq_qual);
    void SetFlattenedPosition(CONTIG_OFFSETS const& contig_offsets,
                              CONTIG_OFFSETS::const_iterator * contig_iter);

    static void SetGlobalFlags(SAM_QNAME_FORMAT _qname_format);
    static SAM_QNAME_FORMAT sam_qname_format;
    static size_t (* parse_fragment_id)(char const* qname);

    SamLine();

    ~SamLine();


    //bool operator<(SamLine const& b) const;
    //bool operator==(SamLine const& b) const;

    SAM_PARSE load(FILE * seqfile);

    void print(FILE * samfile, bool flip_query_strand_flag) const;

    void print_fastq(FILE * fastq_fh) const;

    int alignment_score(char const* tag, int default_if_missing, bool * has_score) const;
    char alignment_space(char default_if_missing, bool * has_space) const;

    //queries the samline for the desired tag code, loading type and value_string if found.
    //returns whether the tag was found
    bool has_tag(char const* tag_code, char & type, char * value_string) const;


    //adds the desired tag. 
    void add_tag(char const* tag_code, char type, char const* value_string);

    size_t aligned_read_length();
    size_t raw_read_length() const { return strlen(this->seq); }
    size_t ones_based_pos() const { return this->pos + 1; }
    size_t zero_based_pos() const { return this->pos; }

    bool paired_in_sequencing() const { return (this->flag & SamFlags::PAIRED_IN_SEQUENCING) != 0; }
    bool mapped_in_proper_pair() const { return (this->flag & SamFlags::MAPPED_IN_PROPER_PAIR) != 0; }
    bool query_unmapped() const { return ((this->flag & SamFlags::QUERY_UNMAPPED)) != 0; }
    bool mate_unmapped() const { return (this->flag & SamFlags::MATE_UNMAPPED) != 0; }

    //strandedness flags are nonzero for the negative strand.
    bool query_on_pos_strand() const { return (this->flag & SamFlags::QUERY_ON_NEG_STRAND) == 0; }
    bool mate_on_pos_strand() const { return (this->flag & SamFlags::MATE_ON_NEG_STRAND) == 0; }

    bool first_read_in_pair() const { return (this->flag & SamFlags::FIRST_READ_IN_PAIR) != 0; }
    bool second_read_in_pair() const { return (this->flag & SamFlags::SECOND_READ_IN_PAIR) != 0; }
    bool alignment_not_primary() const { return (this->flag & SamFlags::ALIGNMENT_NOT_PRIMARY) != 0; }
    bool failed_quality_check() const { return (this->flag & SamFlags::FAILED_QUALITY_CHECK) != 0; }
    bool pcr_or_optical_duplicate() const { return (this->flag & SamFlags::PCR_OR_OPTICAL_DUPLICATE) != 0; }

    char const* mate_ref_name() const;

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


/* void print_sam_line(FILE * sam_fh, */
/*                     char const* qname, int flag,  */
/*                     char const* rname, size_t pos, */
/*                     size_t mapq, char const* cigar, */
/*                     char const* mrnm, size_t mpos, */
/*                     int isize, char const* seq, */
/*                     char const* qual, */
/*                     char const* tag_string, */
/*                     bool output_is_ones_based); */

void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh);

#endif // _SAM_HELPER_H

