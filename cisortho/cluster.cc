#include "diversity_cluster.h"
#include "specific_cluster.h"
#include "hit.h"

using std::cout;

namespace cis {

  //parse clustering criterion
  //initialize ClusterCriterion
//   cluster * LoadCluster(VNT const& j){
//     string kind = jf::gets(j, "kind", "no_clustering");
//     cluster * ret = NULL;
//     if (kind == "diversity") ret = new diversity_cluster(j);
//     else if (kind == "specific") exit(50);
//     return ret;
//   }
  
  void printcluster(string desc, RIT_M_R b, RIT_M_R bi, RIT_M_R ei, RIT_M_R e){
    if (b == e) return;
    //bool found = false;
    int nc = std::distance(b, e);
    if (nc > 15) {
      cout<<"abridged cluster ("<<nc<<" hits), start hit:";
      HIT_CR h = (HIT_CR)**b;
      cout<<" "<<h.start<<"("<<h.name()<<","<<h.score<<","<<h.is_cluster<<")"<<endl;
      return;
    }

    cout<<desc;
    for (RIT_M i = b; i != e; i++){
      HIT_CR h = (HIT_CR)**i;
      if (i == bi || i == ei) cout<<"***";
      cout<<" "<<h.start<<"("<<h.name()<<","<<h.score<<","<<h.is_cluster<<")";
    }
    cout<<endl;
    cout.flush();
    return;
  }

  //maintain a set 'eq' that represents ids that should all be considered
  //part of the same cluster.  This does not allow any discrimination between
  //different cluster types...
  void MergeOverlappingClusters(REGIONS_MULTI &hits){
		
    set<int> eq;
    set<int>::iterator sit;
		
    for (RIT hi = hits.begin(); hi != hits.end(); hi++){
      HIT_R h = (HIT_R)**hi;
      set<int> & t = h.cluster_tags;
      if (t.empty()) {}
      else {
        int tmin = *t.begin();
        if (eq.find(tmin) != eq.end()) eq.insert(t.begin(), t.end());
        else { eq.clear(); eq.insert(tmin); }
        int emin = *eq.begin();
        t.clear();
        t.insert(emin);
      }
    }
  }
	

  //reduce the cluster numbers by mapping them
  int ReduceClusterNumbers(REGIONS_MULTI &hits, int init){

    set<int> eq;
    set<int>::iterator sit;
		
    int cnum = init;
    std::map<int, int> hi2lo;
    
    set<int>::iterator tags_it;

    for (RIT hi = hits.begin(); hi != hits.end(); hi++){
      HIT_R h = (HIT_R)**hi;
      for (tags_it = h.cluster_tags.begin();
           tags_it != h.cluster_tags.end(); ++tags_it){
        int c = *tags_it;
        if (hi2lo.find(c) == hi2lo.end()) hi2lo[c] = cnum++;
      }
    }
    
    set<int> temp_set;
    for (RIT hi = hits.begin(); hi != hits.end(); hi++){
      temp_set.clear();
      HIT_R h = (HIT_R)**hi;
      for (tags_it = h.cluster_tags.begin();
           tags_it != h.cluster_tags.end(); ++tags_it){
        int c = *tags_it;
        temp_set.insert(hi2lo[c]);
        h.cluster_tags = temp_set;
      }
    }

    return cnum;
  }


} // namespace cis
