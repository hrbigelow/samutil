#ifndef _DIVERSITY_CLUSTER_H
#define _DIVERSITY_CLUSTER_H

#include <map>
#include <string>

#include "cluster_base.h"

/*
  specifies a set of minimum hits to be found in each cluster group
  and a minimum number of these minimums to be fulfilled
  for example:
  
  promoter 2
  enhancer 2
  aux_factor 1
  
  with a total number of requirements = 2

  this means that, in order to be considered a cluster, any combination of
  two of the above three rules must be fulfilled.  the promoter rule is fulfilled
  when there are two or more hits that are 'promoters'.  The enhancer
  rule is fulfilled with two or more 'enhancer' hits, etc.
  The variable 'minimum_counts' holds the rules, while the
  minimum_tallies

*/

class diversity_cluster : public cis::cluster {

  typedef std::map<std::string, int> TALLY;
  typedef TALLY::value_type TDAT;
  typedef TALLY::iterator TLIT;

 public:
  diversity_cluster(std::string const& name_, int max_nbor_span_,
                    TALLY const& min_counts_, int min_reqs_);

/*   diversity_cluster(VNT const&); */

  void add_stat(cis::RIT_M_R); //add the statistic of one hit to the cluster statistics
  void drop_stat(cis::RIT_M_R); //drop the statistic of one hit
  bool valid(); //evaluate the statistics collected with add_stat to see if it's a valid cluster
  void clear_stats();
  void assign_clusters();
  void printsafe(cis::RIT_M_R);

 protected:
  TALLY min_counts; //set of individual rules
  TALLY score_tally;
  int min_reqs; //minimum number of requirements to satisfy
  int num_density_reqs; //set to minimum_score_hits.size()

};

#endif // _DIVERSITY_CLUSTER_H
