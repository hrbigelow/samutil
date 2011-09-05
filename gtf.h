#ifndef _GTF_H
#define _GTF_H

#include <cstdlib>
#include <string>
#include <map>
#include <set>

#include "seq_projection.h"

#include "cisortho/memory.h"
#include "cigar_ops.h"
#include "dep/tools.h"

class GTFEntry
{
 public:

    char line[10000];
    char seqname[1000];
    char source[1000];
    char feature[1000];
    size_t start;
    size_t end;
    float score; // 0 to 1000
    char strand; // '+', '-', or '.' (don't care)
    char frame; // '0', '1', '2', or '.' (not coding)
    char attribute_string[1000];

    char gene_id[1000];
    char transcript_id[1000];

    bool is_data_line;

    bool get_next_record(FILE * gtf_fh);
    void print(FILE * gtf_fh);
    void clear_record();

    inline size_t start_bound() const { return this->start - 1; }
    inline size_t end_bound() const { return this->end; }
    inline size_t length() const { return this->end - this->start + 1; }
};



std::map<std::string, std::string>
gtf_to_transcript_gene_map(char const* gtf_file);

std::set<SequenceProjection>
gtf_to_sequence_projection(char const* gtf_file, char const* species);


struct LessSequenceProjectionTarget
{
    bool operator()(SequenceProjection const& a,
                    SequenceProjection const& b) const;
};


struct unique_transcript
{
    std::string contig_name;
    std::string transcript_id;
    char strand;

    unique_transcript(char const* _cn, char const* _tid, char _st);
    unique_transcript();

    bool operator<(unique_transcript const& ut) const;
};


struct TranscriptInfo
{
    float relative_level;
    size_t count;
    size_t length;
    TranscriptInfo(float _rl, size_t _ct, size_t _ln) : 
        relative_level(_rl), count(_ct), length(_ln) { }
    TranscriptInfo() : relative_level(0.0), count(0), length(0) { }
};


typedef std::map<unique_transcript, TranscriptInfo> TRANSCRIPT_EXPR;



#endif // _GTF_H
