#ifndef _SEARCH_H
#define _SEARCH_H


#include <vector>
#include <utility>
#include <iostream>

//#include "dna_types.h"
#include "region.h"
#include "dnacol.h"

namespace cis { 
  class pattern; 
}

class litestream;

struct hit_info;
struct scan_interval;

template <typename T> class nuc_trie;
template <typename T> class index_trie_scan;

/* template <typename TAG> class nuc_trie; */
typedef nuc_trie<hit_info> hit_trie;

typedef std::pair<hit_trie *, scan_interval> LOCUS;
typedef std::vector<LOCUS> LOCI;

LOCI Search(index_trie_scan<hit_trie> &, hit_trie *, 
            litestream &, litestream &, int);

int64_t TotalSize(index_trie_scan<hit_trie> & scantrie,
                    litestream & trie_stream);

int CountHits(LOCI const&);

void FindThresholdScore(LOCI const& loci, 
                        int nhits, std::pair<int, int> * threshold);

cis::REG_MAP makeHits(LOCI const& ivals,
					  cis::pattern *pat,
					  litestream & pos_file,
					  cis::dna_collection const& dnac,
					  cis::dna_collection::const_iterator start_dna,
					  cis::dna_collection::const_iterator end_dna,
					  int max_depth);

#endif //_SEARCH_H
