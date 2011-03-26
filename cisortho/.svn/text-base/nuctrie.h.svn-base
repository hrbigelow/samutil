#ifndef _NUC_TRIE
#define _NUC_TRIE

#include <limits.h>

#include <map>
#include <set>

#include "dna.h"

using std::map;
using std::set;
using cis::NUC;

/* this is the API of a (group of) patterns represented by a trie to
   be used in the 'scan' function.  

*/


//base class for creating a nucleotide trie.  Defines an API for building
//traversing, and collecting labels from the trie.
template <typename TAG>
class nuc_trie {

 public:

  typedef map<NUC, nuc_trie<TAG> *> MAP;
  typedef typename MAP::iterator IT;

  virtual bool leaf() const = 0;
  virtual int score(cis::dnastring const& hit_sequence) const = 0;
  virtual map<NUC, nuc_trie<TAG> *> children(NUC) = 0;
  virtual void collect(set<nuc_trie<TAG> *> &, NUC const*const, int) = 0;
  virtual void reset() = 0;
  virtual int depth() const = 0;
  virtual void apply_to_leaves(void (fcn)(void * accu, NUC *pfx, int sz), 
                               void * accu, int pind = -1);

  virtual TAG const* get_label() const = 0;
  virtual ~nuc_trie() {};

 protected:

};




//the problem is we can only know the base of a given node by its parent
//so, we have to set the prefix buffer before traversing it...
//traverse all nodes of the tree, apply the function to them on the way up
template <typename TAG>
void nuc_trie<TAG>::apply_to_leaves(void (fcn)(void *accu, NUC *pfx, int sz), 
									void * accu, int pind){

  //get all kids consistent with the score thresholds set.
  ++pind;

  static NUC prefix_buffer[1000];

  if (leaf()) {
    fcn(accu, prefix_buffer, pind);
    return;
  }

  MAP kidsA = children(cis::A);
  MAP kidsC = children(cis::C);
  MAP kidsG = children(cis::G);
  MAP kidsT = children(cis::T);

  MAP kids;
  kids.insert(kidsA.begin(), kidsA.end());
  kids.insert(kidsC.begin(), kidsC.end());
  kids.insert(kidsG.begin(), kidsG.end());
  kids.insert(kidsT.begin(), kidsT.end());

  for (IT it = kids.begin(); it != kids.end(); ++it){
    prefix_buffer[pind] = (*it).first; //set the prefix
    (*it).second->apply_to_leaves(fcn, accu, pind);
  }
}


struct hit_info {
  string id;
  cis::dna_strand strand;
  int size;
  int score;
	
  hit_info(string i, cis::dna_strand d, int sz) : id(i), strand(d), size(sz) {}
  hit_info() : id(""), strand(cis::NON_STRANDED), size(0) {}
  friend bool operator<(hit_info const& a, hit_info const& b) {
    return 
      a.id < b.id ||
      (a.id == b.id && 
       a.strand < b.strand);
  }
	
};


const int QUAL_SCORE = INT_MIN + 1;
const int MIN_SCORE = INT_MIN;
const int MAX_SCORE = INT_MAX;


typedef nuc_trie<hit_info> hit_trie;


#endif //_NUC_TRIE
