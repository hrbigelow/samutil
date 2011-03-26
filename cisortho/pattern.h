#ifndef _PATTERN_H
#define _PATTERN_H

#include <string>
#include <vector>
#include "nuctrie.h"
#include "pwm.h"


namespace cis {

class pattern {
  
 public:

  void * ptr;
  void * ptr2;
  hit_info label;
  nuc_trie<hit_info> * t;
  std::string kind;
  std::string cluster_group;
  int max_length;
  int total_hits;
  int min_score;
/*   pattern(VNT const&); //create a pattern from a parsed JSON variant type */
  pattern(void * ptr, void * ptr2, hit_info label, nuc_trie<hit_info> * ntrie,
          std::string const& kind, std::string const& cg, int max_length,
          int total_hits, int min_score);

  pattern(pattern const&);
/*   void init(VNT const&); */
  ~pattern();
  void init_matrix(IMATRIX const& matrix);
  void init_collection(std::map<std::string, cis::dnastring> const& 
                       consensus_pats);
  
};

 pattern * MakeRCPattern(pattern const* pat);
   
   /*  pattern * JSONtoPattern(VNT const& s); */
/*  pattern * JSONtoPatternRC(VNT const& s); */

} // namespace cis

#endif // _PARSE_H
