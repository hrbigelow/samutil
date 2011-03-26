#include "diversity_cluster.h"
#include "region.h"
#include "hit.h"

using namespace cis;

/* Simply require a minimum number of clusters with certain scores... */
// diversity_cluster::diversity_cluster(char const* diversity_constraint, 
// 									 int nbor_span) :

diversity_cluster::diversity_cluster(std::string const& name_, int max_nbor_span_,
                                     diversity_cluster::TALLY const& min_counts_,
                                     int min_reqs_) :
  cluster(name_, max_nbor_span_),
  min_counts(min_counts_), min_reqs(min_reqs_) {
  _cluster_bitmask = 1<<cluster_count++;
}



// diversity_cluster::diversity_cluster(VNT const& cl) {

// 	name = jf::gets(cl,"name");
// 	max_neighbor_span = jf::geti(cl,"window");
// 	min_counts = ParseHash(jf::getv(cl,"reqs"), &jf::Kfcn<std::string>, &jf::exti);
// 	min_reqs = jf::geti(cl,"min_reqs");
// 	_cluster_bitmask = 1<<cluster_count++;
// }



void diversity_cluster::add_stat(RIT_M_R hi){

	HIT_PC h = (HIT_PC) *hi;
	//std::string & cg = *h->cluster_group;
	++score_tally[*h->cluster_group];
// 	if (score_tally.find(cg) != score_tally.end())
// 		score_tally[cg]++;
// 	int id = h->cluster_group;
// 	for (int a = 0; a < num_density_reqs; a++) score_tally[a][id]++;
}


void diversity_cluster::drop_stat(RIT_M_R hi){

	HIT_PC h = (HIT_PC) *hi;
	--score_tally[*h->cluster_group];
// 	int id = h->cluster_name.uid();
// 	for (int a = 0; a < num_density_reqs; a++) score_tally[a][id]--;
}



bool diversity_cluster::valid(){

	int nqual = 0;
    TALLY::iterator tit;
    for (tit = min_counts.begin(); tit != min_counts.end(); ++tit)
      if (score_tally[(*tit).first] >= (*tit).second) ++nqual;

	return nqual >= min_reqs;

	//c is the constraint
// 	for (int c = 0; c < num_density_reqs; c++){
// 		int nqual = num_algos;
// 		for (int t = 0; t < num_algos; t++)
// 			if (score_tally[c][t] < minimum_counts[c][t]) nqual--;

// 		if (nqual < minimum_tallies[c]) {
// 			qual = false;
// 			break;
// 		}
// 	}

// 	return qual;
}





//uses the tags method (multiple tags per hit, allowing many-to-one membership)
//does this work across several DNAs?  Doesn't do that...just checks
//the start position of the hit, not the DNA...
void diversity_cluster::assign_clusters(){

	//scan the window incrementally
	RIT_M e = hits->end();

	clear_stats();

	RIT_M lf = hits->begin();
	RIT_M rt_cl = lf;
	RIT_M rt_cl_old = lf;
	RIT_M lf_cl = lf;

	while (rt_cl != e){

		rt_cl = bound(lf_cl, rt_cl_old, e, max_neighbor_span);
		for (RIT_M c = rt_cl_old; c != rt_cl; ++c) add_stat(c);
		if (valid()) {
			add_tag(lf_cl, rt_cl);
			set_cluster_flag(lf_cl, rt_cl);
		}

		if (rt_cl == e) break;
		drop_stat(lf_cl);
		++lf_cl;
		rt_cl_old = rt_cl;
	}

}


void diversity_cluster::clear_stats(){
	score_tally.clear();
}
