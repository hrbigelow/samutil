#ifndef _NESTED_H
#define _NESTED_H

#include <utility>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>

struct ContigCoord
{
    std::string contig;
    size_t coord;
    ContigCoord(char const* _contig, size_t _coord) : 
        contig(_contig), coord(_coord) { }

    bool operator<(ContigCoord const& b) const
    {
        int contig_cmp = strcmp(this->contig.c_str(), b.contig.c_str());
        return contig_cmp < 0
            || (contig_cmp == 0
                && this->coord < b.coord);
    }
};

class IntervalTree;

typedef std::map<ContigCoord, IntervalTree *> BY_START;

  
class IntervalTree {
    
 public:
    typedef std::pair<size_t, size_t> INTERVAL;

	IntervalTree(char const* contig, size_t start, size_t end, void * payload);
    IntervalTree();
	~IntervalTree();
    bool has_interval;
    char * contig; //owns the contig.  the top-level node does new/delete for it.
    size_t start;
    size_t end;
    void * payload; //does not own

    BY_START * children; //will be the owner of the subtrees

	std::pair<BY_START::iterator, bool> 
        insert_rec(char const* contig, size_t start, size_t end, void * payload); 

	std::pair<BY_START::iterator, bool> 
        insert(char const* contig, size_t start, size_t end, void * payload); 

    std::vector<IntervalTree const*> 
        overlap(char const* contig, size_t start, size_t end);

 private:
	IntervalTree(IntervalTree const&);
	IntervalTree & operator=(IntervalTree const&);
	IntervalTree(IntervalTree &);
    IntervalTree & operator=(IntervalTree &);
};

typedef std::vector<IntervalTree::INTERVAL> INTERVALS;


std::pair<BY_START::const_iterator, BY_START::const_iterator> 
RangeQuery(BY_START const*, char const* contig, size_t start, size_t end);

bool node_contains(IntervalTree const*, IntervalTree const*);
bool well_ordered(IntervalTree const*, IntervalTree const*);

bool check_tree(IntervalTree const* root);
void print_tree(IntervalTree const* node, int indent_level);

template<typename Payload> void free_payload(IntervalTree * root)
{
    if (root->payload != NULL)
    {
        delete static_cast<Payload *>(root->payload);
        root->payload = NULL;
    }
    if (root->children == NULL)
    {
        return;
    }
    else
    {
        for (BY_START::iterator child = root->children->begin();
             child != root->children->end(); ++child)
        {
            free_payload<Payload>((*child).second);
        }
    }
}

#endif //_NESTED_H
