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


enum SAM_ORDER
    {
        SAM_RID,
        SAM_POSITION_RID,
        SAM_RID_POSITION,
        SAM_POSITION,
        SAM_FID_POSITION
    };


class SamLine
{
public:
    static bool numeric_start_fragment_ids; 
    // set this to true if we're in a training situation with
    // numerically named qnames.

    size_t bytes_in_line;
    char * line;
    SAM_PARSE parse_flag;
    char * qname;
    size_t qid;
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

    SamLine(SAM_PARSE _parse_flag,
            char const* _qname, int _flag, 
            char const* _rname, size_t _pos,
            size_t _mapq, char const* cigar,
            char const* _mrnm, size_t _mpos,
            int _isize, char const* _seq,
            char const* _qual,
            char const* _tag_string);

    SamLine(FILE * sam_fh, bool sam_is_ones_based, bool allow_absent_seq_qual);
    SamLine(char const* samline_string, bool sam_is_ones_based, bool allow_absent_seq_qual);

    void Init(char const* samline_string, bool sam_is_ones_based, 
              bool allow_absent_seq_qual);

    SamLine(SamLine const&)
    {
        assert(false);
    }
    SamLine();

    ~SamLine();

    //uses fields [rname, pos, strand, cigar]
    bool less_position(SamLine const& b) const;
    bool equal_position(SamLine const& b) const;

    //uses fields [rname, first_read_in_pair()]
    bool less_rid(SamLine const& b) const;
    bool equal_rid(SamLine const& b) const;

    //sort order [qname, pair, rname, pos, query_strand, cigar, mrnm, mpos]
    bool less_rid_position(SamLine const& b) const;
    bool equal_rid_position(SamLine const& b) const;

    // sort order [rname, pos, query_strand, cigar, mrnm, mpos, qname, pair]
    bool less_position_rid(SamLine const& b) const;
    bool equal_position_rid(SamLine const& b) const;

    // sort order [rname]
    bool less_fid(SamLine const& b) const;
    bool equal_fid(SamLine const& b) const;

    // sort order [rname, pos, strand, cigar]
    bool less_fid_position(SamLine const& b) const;
    bool equal_fid_position(SamLine const& b) const;

    //bool operator<(SamLine const& b) const;
    //bool operator==(SamLine const& b) const;

    SAM_PARSE load(FILE * seqfile, bool source_is_ones_based_pos);

    void print(FILE * samfile, bool flip_query_strand_flag) const;

    void print_fastq(FILE * fastq_fh) const;

    int alignment_score(char const* tag, int default_if_missing, bool * has_score) const;
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


struct less_char_ptr
{
    bool operator()(char const* a, char const* b) const
    {
        return strcmp(a, b) < 0;
    }
};

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

//typedef std::unordered_map<char const*, size_t, std::hash<char const*>, eqstr> CONTIG_OFFSETS;
typedef std::map<char const*, size_t, less_char_ptr> CONTIG_OFFSETS;
typedef CONTIG_OFFSETS::const_iterator OFFSETS_ITER;


std::map<std::string, size_t> ContigLengths(FILE * sam_fh);

CONTIG_OFFSETS ContigOffsets(std::map<std::string, size_t> const& contig_lengths);

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


void print_sam_line(FILE * sam_fh,
                    char const* qname, int flag, 
                    char const* rname, size_t pos,
                    size_t mapq, char const* cigar,
                    char const* mrnm, size_t mpos,
                    int isize, char const* seq,
                    char const* qual,
                    char const* tag_string,
                    bool output_is_ones_based);

size_t flattened_position(char const* contig, size_t position, size_t flag,
                          CONTIG_OFFSETS const& contig_offsets,
                          CONTIG_OFFSETS::const_iterator * contig_iter);

size_t samline_position_align(char const* samline, CONTIG_OFFSETS const& contig_offsets);

size_t samline_read_id_flag(char const* samline, CONTIG_OFFSETS const& /* unused */);

size_t index_sam_rid_position(char const* samline, CONTIG_OFFSETS const& contig_offsets);

size_t samline_position_min_align_guide(char const* samline, 
                                        CONTIG_OFFSETS const& contig_offsets);

void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh);

#endif // _SAM_HELPER_H

