#ifndef _HIT_H
#define _HIT_H

#include <limits.h>

#include <iostream>
#include <set>
#include <utility>
#include <string>

#include "dna.h"
#include "region.h"
#include "memory.h"
#include "dna_types.h"

using std::set;
using std::string;

namespace cis {

  typedef std::map<cis::dnastring, unsigned int> SEQMAP;
  typedef std::pair<cis::dnastring, unsigned int> SEQDAT;
  typedef SEQMAP::iterator SEQ_IT;
  
  enum hit_status { ORPHAN, NON_ORPHAN, NOT_KNOWN } ;
  
  const unsigned int UNASSIGNED_CLUSTER = INT_MAX;
  
  //class for storing and retrieving a string value
//should maintain a single static set of stored data
//instances of store will 
 class hit;


typedef hit & HIT_R;
typedef hit const& HIT_CR;
typedef hit * HIT_P;
typedef hit const* HIT_PC;
typedef hit const*const HIT_PCP;


class hit : public region {

	static unsigned int seq_id_ctr;

 public:

	//unsigned int index;
	memory::cstore<cis::dnastring> seq;
	//cis::dnastring const& seq;
	int const score;
	//REG_P base; //identifies the containing region...
	//scoring const* scoring_alg;
	memory::cstore<string> name;
	//store<string> cluster_name;
	string * cluster_group; //the equivalence class that this hit belongs to for clustering...
	set<int> cluster_tags; //this must discern clusters of different types...
	int is_cluster; // bitset of clusters this hit belongs to

	//string * sub_name; //to be used for TRIE-based scoring, where the algorithm
	//name is not sufficient

	friend ostream& operator<<(ostream&, HIT_CR);
	friend istream& operator>>(istream&,hit&);
	friend bool operator==(HIT_CR, HIT_CR);
	friend bool operator<(HIT_CR, HIT_CR);
	hit(dna_t const&, cis::dnastring const& sequence, 
        uint64_t start, uint64_t end, 
		dna_strand const& strand, int score,
		int id = -1, string name = "", 
        string * cluster_group = NULL);
	~hit();

	string printbrief() const;
	//string name() const;

	static unsigned int hitcounter;
	static void SetCounter(unsigned int i=0){ hitcounter = i; }
	//static cis::dnastring const& store_seq(cis::dnastring const&);

};


 void PrintHitClusters(REGIONS_MULTI const& hits, FILE *outstream);

 void PrintHits(REGIONS_MULTI const& hits, dna_t const& dna, FILE *outstream);



struct less_hit_score_ptr { bool operator()(HIT_PC, HIT_PC) const; };

typedef std::set<HIT_P, less_hit_score_ptr> HITS_SORT;

typedef HITS_SORT::iterator SIT;

std::pair<REG_PC, REG_PC> FindBoundingExons(REGIONS_MULTI const&, HIT_PC);

} //namespace cis

#endif
