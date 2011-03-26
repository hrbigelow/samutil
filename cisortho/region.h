#ifndef _REGION
#define _REGION

#include <map>
#include <utility>
#include <string>
#include <sstream>

#include "dna.h"
#include "enum.h"
#include "region_types.h"

const int GENERIC_REGION_ID = -2;


//these types are relative to a single gene
//have a class that simulates an O'Caml polymorphic type...

namespace cis {
  class dna_t;
  typedef std::set<dna_t const*, less_dna_ptr> DNAS;
  class region_tree;
}


namespace cis {

  std::string PrintEnum(int d);

  class exon; 
  typedef exon const* EXON_PC;
  typedef exon const& EXON_RC;
  typedef exon * EXON_P;

  class region {

  public:
    cis::dna_t const& dna;
    int id;
    int group_id;
    int64_t start, end;
    dna_strand strand;
    std::string name;
    int group_relation_;
    void const* payload;
    inline int64_t length() const { return end - start; } //assumes end >= start
    region(cis::dna_t const& dna, std::string name, int64_t start, 
           int64_t end, dna_strand strand, int id, int group,
           int group_relation, void const* payload);
/*     region(cis::dna_t const&, std::string const& = "", int64_t = 0, int64_t = 0, */
/*            dna_strand = NON_STRANDED, int = -1, int = 0, int gr = region::NONE); */
    bool contains(REG_CR) const;
    bool overlaps(REG_CR) const;
    bool same_strand(EXON_PC) const;
    cis_enum::rel_position relative_position(EXON_PC) const;
    int offset(EXON_PC) const;
/*     cis_enum::syn_type synteny_type() const; */
    friend bool operator<(region const&, region const&);
    friend bool operator==(region const&, region const&);
    static int region_id;

    static const int NONE;
    static const int LEFT_MOST;
    static const int RIGHT_MOST;
  };

  struct less_region_ptr { bool operator()(REG_P, REG_P) const; };

  struct less_region_pos {
    bool operator()(REG_PC a, REG_PC b) const { 
      return 
        a->start < b->start ||
        (a->start == b->start &&
         (a->end < b->end ||
          (a->end == b->end &&
           a->id < b->id)));
    }
  };

/*   typedef multi_index_container< */
/*     REG_PC, */
/*     indexed_by< */
/*     ordered_unique<tag<start_tag>, identity<REG_PC>, less_region_pos>, */
/*     ordered_unique<tag<index_tag>, BOOST_MULTI_INDEX_MEMBER(region, const int, id)> */
/*     > */
/*     > REGIONS_MULTI; */

  typedef std::set<region const*, less_region_pos> REGIONS_MULTI;

  typedef REGIONS_MULTI::iterator RIT;
  typedef REGIONS_MULTI::iterator RIT_M;
  typedef REGIONS_MULTI::iterator const& RIT_M_R;

  typedef REGIONS_MULTI::reverse_iterator RIT_M_REV;
  typedef REGIONS_MULTI::reverse_iterator const& RIT_M_REV_R;

  typedef REGIONS_MULTI::reverse_iterator RIT_REV;

  typedef std::map<cis::dna_t const*, REGIONS_MULTI> REG_MAP;
  typedef std::pair<cis::dna_t const*, REGIONS_MULTI> REG_PAIR;


  int64_t RegionDistance(REG_CR r1, REG_CR r2);
  int64_t RegionOverlap(REG_CR r1, REG_CR r2);


  struct less_region_start {
    bool operator()(REG_PC a, REG_PC b) const { 
      return 
        a->start < b->start ||
        (a->start == b->start &&
         (a->end < b->end ||
          (a->end == b->end &&
           a->id < b->id)));
    }
  };


  struct less_region_end {
    bool operator()(REG_PC a, REG_PC b) const { 
      return 
        a->end < b->end ||
        (a->end == b->end &&
         (a->start < b->start ||
          (a->start == b->start &&
           a->id < b->id)));
    }
  };



/*   REGIONS_MULTI overlapping(REGIONS_MULTI const& reg, REG_PC r); */
/*   std::string lower_case_regions(std::string const& seq, REG_PC r,  */
/*                             cis::region_tree * regs); */

  //now, instead of getting the left bound this way, just return the
  //reverse iterator and take its base if you want the left bound...


  template<typename FI>
    FI bound(FI const& start,
             FI const& hint,
             FI const& ub, 
             int const& max_dist){
	
    FI c = hint;
    int lmid = (*start)->start;
	
    if (c == ub || (abs((*c)->start - lmid)) > max_dist) return c;
    else {
      do c++;
      while (c != ub && (abs((*c)->start - lmid)) <= max_dist);
      return c;
    }
  }


  template<typename FI>
    FI bound(FI const& start,
             FI const& ub, 
             int const& max_dist){
	
    FI c = start;
    int lmid = (*start)->start;

    if (c == ub || (abs((*c)->start - lmid)) > max_dist) return c;
    else {
      do c++;
      while (c != ub && (abs((*c)->start - lmid)) <= max_dist);
      return c;
    }
  }


  RIT_M_R min_iter(REGIONS_MULTI * c, RIT_M_R a, RIT_M_R b);

  REG_MAP RDBToRegions(DNAS& dnas, char const* rdbfile, bool with_group_field);

  REGIONS_MULTI MergeRegions(REG_MAP r);
  REG_MAP SplitRegions(REGIONS_MULTI const& regions);

  REG_P MakeRegion(cis::dna_t const*, int64_t, int64_t, std::string const&,
                   std::string const&, int, int, int, int = -1, int = 0);

  void SetGroupRelation(RIT const& begin, RIT const& end);
  void PrintDNARegion(REG_PC dna_region, FILE * outfile);
  void PrintDNARegionSimple(REG_PC dna_region, int max_sequence_length, 
                            FILE * outfile);
  
} // namespace cis

#endif //_REGION
