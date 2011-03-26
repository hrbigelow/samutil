#ifndef _CLUSTER_H
#define _CLUSTER_H

#include <string>

#include "region.h"
#include "cluster_base.h"

namespace cis { class cluster; }

namespace cis {

/*   cluster * LoadCluster(VNT const&); */
  void clusterstring(RIT_M_R, RIT_M_R lf, RIT_M_R lfi, RIT_M_R rti, RIT_M_R rt);
  void printcluster(std::string desc, RIT_M_R, RIT_M_R, RIT_M_R, RIT_M_R);
  void MergeOverlappingClusters(REGIONS_MULTI &hits);
  int ReduceClusterNumbers(REGIONS_MULTI &hits, int init);

}

#endif // _CLUSTER_H
