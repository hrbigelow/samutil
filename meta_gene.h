#ifndef _META_GENE_H
#define _META_GENE_H


#include "nclist.h"
#include "cigar_ops.h"

#include <map>
#include <string>

class GTFEntry;

struct FeatureJumps
{
    std::map<size_t, std::pair<int, std::set<std::string> > > data;
    std::map<size_t, std::pair<int, std::set<std::string> > >::const_iterator iter;

    FeatureJumps()
    {
        this->iter  = data.begin();
    }

    void insert(size_t boundary_position, bool is_start, char const* feature_name)
    {
        if (this->data.find(boundary_position) == this->data.end())
        {
            data[boundary_position] = std::pair<int, std::set<std::string> >();
        }
        std::pair<int, std::set<std::string> > & jump = data[boundary_position];
        jump.first += is_start ? 1 : -1;
        jump.second.insert(std::string(feature_name));
    }
    size_t boundary() const
    {
        return this->iter->first;
    }

    int height() const
    {
        return this->iter->second.first;
    }
    std::set<std::string> const& features() const
    {
        return this->iter->second.second;
    }
    void init() { this->iter = data.begin(); }
    void next() { ++this->iter; }
    bool is_last() { return this->iter == data.end(); }
};


struct unique_gene_description
{
    std::string contig;
    std::string gene;
    char strand;
    unique_gene_description(char const* _contig, char const* _gene,
                            char _strand) :
        contig(_contig), gene(_gene), strand(_strand) { }

    unique_gene_description() : contig(""), gene(""), strand('+') { }
    
    bool operator<(unique_gene_description const& u) const;
};



struct BoundsInfo
{
    size_t start_boundary;
    size_t end_boundary;
    size_t start_cigar_index;
    BoundsInfo(size_t _sb, size_t _eb, size_t _sci) : 
        start_boundary(_sb), end_boundary(_eb), start_cigar_index(_sci) { }
    BoundsInfo() :
        start_boundary(0), end_boundary(0), start_cigar_index(0) { }
};



typedef std::map<unique_gene_description, BoundsInfo> GENE_INDEX;
typedef std::set<std::string> CHAR_SET;
typedef std::map<std::string, CHAR_SET> GENE_TRANSCRIPT_MAP;
typedef std::vector<GENE_INDEX::iterator> GENE_ITERS;

class MetaGene
{

 public:
    IntervalTree * gene_extent_tree;
    IntervalTree * meta_exon_tree;

    GENE_INDEX gene_bounds;

    std::map<std::string, Cigar::CIGAR_VEC> meta_exon;
    std::map<std::string, std::multiset<std::pair<size_t, size_t> > > meta_exon_offsets;

    bool bounds_initialized;
    GENE_TRANSCRIPT_MAP transcripts;

    MetaGene();
    ~MetaGene();
    void Initialize(FILE * gtf_fh, size_t pseudo_intron_length);
    size_t transcript_number(GTFEntry const& gtf_entry) const;
    
};

#endif // _META_GENE_H
