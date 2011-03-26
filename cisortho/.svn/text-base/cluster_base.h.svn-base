#ifndef _CLUSTER_BASE_H
#define _CLUSTER_BASE_H

#include <string>

#include "region.h"

namespace cis {

  class cluster {
  public:
    std::string name;
    enum cluster_kind { 
      MIN_DIVERSITY, DIVERSITY, SPECIFIC, NO_CLUSTERING 
    } kind;
		
    static int cluster_count;
    REGIONS_MULTI * hits;

    bool single_tag_mode;

    int max_neighbor_span;
    int _cluster_bitmask;

    cluster(std::string n = "", int span = 0) : name(n), max_neighbor_span(span) {}

    void add_tag(cis::REGIONS_MULTI::iterator const& start, 
                 cis::REGIONS_MULTI::iterator const& end); //add a new tag to the range [start, end)
    void clear_tag(cis::REGIONS_MULTI::iterator const& start, 
                   cis::REGIONS_MULTI::iterator const& end, int tag); //delete tag from all hits in the range [start, end)
    void set_cluster_flag(cis::REGIONS_MULTI::iterator const&, 
                          cis::REGIONS_MULTI::iterator const&);
    void clear_flags(cis::REGIONS_MULTI::iterator const&, 
                     cis::REGIONS_MULTI::iterator const&);
    void clear_tags(cis::REGIONS_MULTI::iterator const&, 
                    cis::REGIONS_MULTI::iterator const&);

    virtual void assign_clusters() = 0;
    bool complete_cluster(cis::REGIONS_MULTI::iterator const&, 
                          cis::dna_t const&, cis::dna_t const&);
    bool complete_span(cis::REGIONS_MULTI::iterator const&, 
                       cis::REGIONS_MULTI::iterator const&);
    std::string printsafe(cis::REGIONS_MULTI::iterator const&);

    void set_hits_container(REGIONS_MULTI *);

    virtual ~cluster() {}

  };

} //namespace cis

#endif //_CLUSTER_BASE_H
