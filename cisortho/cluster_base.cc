#include <cmath>
#include <set>

#include "cluster_base.h"
#include "hit.h"

using std::cout;

//matches zic+{-2,30}ebox

namespace cis {

	unsigned int UniqueID(){
		static unsigned int counter = 0;
		return counter++;
	}

	int cluster::cluster_count = 0;

	//find the hit to the left of <right> such that
	//[<left>, <right>] is the largest range contained in max_neighbor span.  (notice the [] range rather than [) for 
	//right_bound.
	//will return begin() if the

	void cluster::set_cluster_flag(RIT_M_R start, RIT_M_R end){

		for (RIT_M p = start; p != end; p++)
			((HIT_P)*p)->is_cluster |= _cluster_bitmask;
	}
	

	//clear cluster tags in [lf, rt) depending on the prevalent mode
	void cluster::clear_tags(RIT_M_R lf, RIT_M_R rt){
		for (RIT_M p = lf; p != rt; p++) (HIT_P(*p))->cluster_tags.clear();
	}

	void cluster::clear_flags(RIT_M_R lf, RIT_M_R rt){
		for (RIT_M p = lf; p != rt; p++)
			(HIT_P(*p))->is_cluster ^= _cluster_bitmask;
	}
	

	//adds a new tag to all hits in the range [start, end)
	void cluster::add_tag(RIT_M_R start, RIT_M_R end){
	
		if (start == end){
			cout<<"in add tag: start and end are equal!"<<endl;
		}

		int tag = UniqueID();
		//cout<<"adding "<<tag<<" to hit_id range ["<<(*start)->index<<","<<(*end)->index<<")"<<endl;
		for (RIT_M h = start; h != end; h++) (HIT_P(*h))->cluster_tags.insert(tag);
	}



	void cluster::clear_tag(RIT_M_R start, RIT_M_R end, int tag){
		for (RIT_M h = start; h != end; h++) (HIT_P(*h))->cluster_tags.erase(tag);
	}




	/************************************************************************
assumptions for cluster assignment:
the range given is guaranteed to be a bounded set.  That is,
no definition of cluster would otherwise cross the boundary defined by the given input.
second:  the function is allowed to define an arbitrary number of clusters by inserting
values into the cluster tags taken from a unique ID generator...
	************************************************************************/




	/*
	  Updates the cluster status of [lfi,rti).  If lfi is part of a cluster,
	  it means that cluster may reside in [left_bound(lfi), lfi+1).  likewise,
	  if rti-1 is part of a cluster, that cluster may reside in [rti-1, right_bound(rti-1))
	  this larger extent is that which must be scanned to re-evaluate the cluster
	  status of [lfi,rti).

	  If re_eval is false, treats the range [lfi,rti) such that it assumes a barrier of
	  independence at lfi and rti.  That is, there is no cluster
	*/
	string cluster::printsafe(RIT_M_R r){

		std::ostringstream os;
		if (r == hits->end()) os<<"end";
		else os<<(*r)->start<<"(id "<<(*r)->id<<")";
		return os.str();
	}


	void cluster::set_hits_container(REGIONS_MULTI *h) { hits = h; }




} // namespace cis
