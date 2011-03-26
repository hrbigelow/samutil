#include "nuctrie.h"

#include <string>
#include "dna.h"
#include "hit.h"

using cis::NUC;
using cis::dna_strand;
using cis::dnastring;
using std::string;
using cis::hit_info;

/*
  for searching a prefix trie with a consensus pattern, allowing
  up to m mismatches, which will be recorded as score = m - # mismatches.
  The node must simulate the behavior of a tree representing a collection
  of patterns and their scores.  This means that each node must have a constant
  score attached to it.  
  
  Should the label contain all hit-specific information, or should the path
  be stored in the 
 */

class consensus_trie : public nuc_trie<hit_info> {
	
	consensus_trie * next;
	NUC _base;
	NUC _childbits;
	hit_info * _label;
	int mismatch;
	int max_mismatch;

	consensus_trie(int m, int mm, NUC n);
	~consensus_trie();
	NUC childbits();
	int score();
	consensus_trie * child(NUC);
	void insert(NUC const*const pat, int, hit_info const & label);
	void insert(dnastring, hit_info const & label);
	void collect(set<nuc_trie<hit_info> *> &, NUC const*const, int);
	void collectm(set<nuc_trie<hit_info> *> &, NUC const*const, int);

};
