#ifndef _NESTED_H
#define _NESTED_H

#include <utility>
#include <map>
#include <set>

#include "region.h"

namespace cis {

  typedef std::map<int64_t, region_tree *> BY_START;
  typedef std::map<int64_t, region_tree *> BY_END;
  
  class dna_t;
  
  class region_tree {

    
  public:

	region_tree(region const* r);
	~region_tree();
	region const* interval;
    BY_START * children; //will be the owner of the subtrees
/*     BY_END * children_by_end; */
	void insert_rec(region const*); //creates a new node and inserts it appropriately into the tree structure...
	void insert(region const*); //creates a new node and inserts it appropriately into the tree structure...

  private:
	region_tree(region_tree const&);
	region_tree & operator=(region_tree const&);
	region_tree(region_tree &);
	region_tree & operator=(region_tree &);
	region_tree();
  };

  typedef std::map<dna_t const*, region_tree *> TREE_MAP;


  REGIONS_MULTI IntervalOverlap(region_tree const&, region const&);

  std::pair<BY_START::const_iterator, BY_START::const_iterator> 
    RangeQuery(BY_START const*, region const&);

  std::pair<region const*, region const*> 
    NonOverlappingNeighbors(BY_START const*, BY_END const*, region const& q);

  region_tree * BuildTree(REGIONS_MULTI const&);
  TREE_MAP BuildTrees(REGIONS_MULTI const& regions);

  
  std::set<std::pair<region const*, region const*> > 
    GetAllNeighbors(REGIONS_MULTI const&, region_tree const&, int);
  

  //a safe and efficient test of less
  //template <typename ITER, typename access> 
  bool less_iterator(BY_START::iterator & a, BY_START::iterator & b, 
                     BY_START::iterator & end);

  void print_tree(region_tree const* node, int indent_level);


  /*
    int node_index (node const& n){
	return n.interval.start;
    }


    struct less_node_start {
	bool operator()(node const& a, node const& b){
    return node_index(a) < node_index(b);
	}
    };


    typedef std::set<node, less_node_start> NSET;


    //what should be the interface for the nested class?
    //besides insert and erase, size, empty, what other things?
    //we want to be able to do a range query.
    class nested : public NSET {};
	



    //ordering consistent with the forward depth first search traversal
    //of a nested containment tree structure
    struct forward_depth_first {
	bool operator()(node const& an, node const& bn){
    INTERVAL &a = an.interval;
    INTERVAL &b = bn.interval;
    return 
    a.start < b.start ||
    (a.start == b.start &&
    (a.end > b.end ||
    (a.end == b.end &&
    (a.id < b.id))));
	}
    };


    //ordering consistent with the depth-first traversal of the same
    //tree, but 
    struct reverse_depth_first {
	bool operator()(node const& an, node const& bn){
    INTERVAL &a = an.interval;
    INTERVAL &b = bn.interval;
    return 
    a.end < b.end ||
    (a.end == b.end &&
    (a.start < b.start ||
    (a.start == b.start &&
    (a.id > b.id))));
	}
    };
  */

} // namespace cis

#endif //_NESTED_H

