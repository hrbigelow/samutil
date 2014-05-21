#ifndef _SAM_INDEX_H
#define _SAM_INDEX_H

#include <cstddef>
#include <cstring>
#include <unordered_map>
#include <map>

class SequenceProjection;

// here, the following nicknames are used:
// FID: "fragment ID"
// ALN: "alignment position"
// MAGALN: "minimum between alignment and guide alignment"
// PROJALN: "projected alignment position"
enum SAM_INDEX_TYPE
    {
        SAM_INDEX_FID = 0,
        SAM_INDEX_FID_ALN,
        SAM_INDEX_FID_MAGALN,
        SAM_INDEX_FID_PROJALN,
        SAM_INDEX_ALN_FID,
        SAM_INDEX_MAGALN_FID,
        SAM_INDEX_PROJALN_FID,
        SAM_INDEX_UNDEFINED
    };

SAM_INDEX_TYPE sam_index_type_from_string(char const* index_type);


enum SAM_QNAME_FORMAT
    {
        SAM_NUMERIC,
        SAM_ILLUMINA,
        SAM_CASAVA18,
        SAM_FORMAT_UNDEFINED
    };


// F = flowcell, L = lane, T = tile, X, Y, P = alignment position
struct flowcell_position_idx_t {
    size_t F : 16;
    size_t L : 8;
    size_t T : 16;
    size_t X : 24;
    size_t Y : 24;
    size_t P : 40;
};

struct position_flowcell_idx_t {
    size_t P : 40;
    size_t F : 16;
    size_t L : 8;
    size_t T : 16;
    size_t X : 24;
    size_t Y : 24;
};


union idx_t {
    flowcell_position_idx_t fp;
    position_flowcell_idx_t pf;
    size_t raw[2];
};

struct less_str
{
    bool operator()(char *a, char *b)
    {
        return strcmp(a, b) < 0;
    }
};

typedef std::map<char *, unsigned int, less_str> index_dict_t;


// new sam_index
struct sam_index
{

    idx_t key;
    size_t start_offset : 48; // offset from a raw buffer holding contents of SAM or BAM records
    size_t line_length : 16; // length of this line in the raw buffer.

    sam_index(size_t * _index, size_t _so, int _ll);
    sam_index();
};


bool samidx_less_offset(sam_index const& a, sam_index const& b);
bool samidx_less_key(sam_index const& a, sam_index const& b);
bool samidx_equal_key(sam_index const& a, sam_index const& b);


struct to_integer : public std::hash<std::string>
{
    size_t operator()(char const* k) const;
};

struct eqstr
{
    bool operator()(const char* s1, const char* s2) const;
};


typedef std::unordered_map<char const*, size_t, to_integer, eqstr> HMAP;

// defines attrubtes of genomic contigs, and optionally projections of a transcriptome
struct contig_dict
{
    // fast lookup of integers assigned to each contig (genomic or transcriptome)
    HMAP name_map;
    size_t * length;
    size_t * offset;
    char * space;

    HMAP projected_name_map;
    SequenceProjection * projections;

    contig_dict();
    ~contig_dict();
};




// initializes lengths, offsets, and name_map fields of dict
void init_contig_length_offset(char const* sam_header, contig_dict * dict);


// initializes projected_name_map and projections fields in dict
void init_sequence_projection(char const* gtf_file, contig_dict * dict);


// initialize the space field of dict
// void init_space();


// set the index fields of idx for this sam_line according to the chosen type and format
void set_sam_index(char const* samline, SAM_INDEX_TYPE itype, SAM_QNAME_FORMAT qfmt, 
                   contig_dict const* dict, 
                   index_dict_t *flowcell_dict,
                   sam_index * idx);


// helper function
size_t samidx_flattened_position(char const* contig, size_t position, contig_dict const* dict);


// update this index's flowcell id in place
void sam_update_flowcell_id(const unsigned int * remap, SAM_INDEX_TYPE itype, sam_index * idx);


sam_index samidx_make_max();

SAM_QNAME_FORMAT qname_format(char const* sam_dataline);


#endif // _SAM_INDEX_H
