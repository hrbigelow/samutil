#include "hit.h"
#include "string_tools.h"
#include <sstream>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

using namespace std;
using namespace cis;

namespace cis
{

	unsigned int hit::hitcounter;
	
	//can we assume the source region has a non-null gene?
	//what we really want is a "region with" ocaml type construct
	hit::hit(dna_t const& dna, dnastring const& s, uint64_t start, uint64_t end, 
			 dna_strand const& strand, 
			 int score, int id, string n, string * cg):
		region(dna, "", start, end, strand, id, -1, region::NONE, NULL), 
		score(score), cluster_group(cg), is_cluster(0) {
		if ((int)s.size() < this->length()){
			cerr<<"Warning: hitstring is shorter than length of hit."<<endl;
		}
		seq(dnastring(s, 0, this->length()));

		name(n);
	}
	
	
	hit::~hit(){
// 		if (sub_name != NULL) { 
// 			delete sub_name; 
// 			sub_name = NULL;
// 		}
	}
	

//calculates the offset of the hit with the gene:
//if hit is outside of gene, the offset is absolute distance
//to nearest exon


//operators for hit
bool operator<(HIT_CR h1,HIT_CR h2){
	return 
		h1.score < h2.score ||
 		(h1.score == h2.score &&
		 ((region&)h1 < (region&)h2));
}

bool operator==(HIT_CR h1,HIT_CR h2){
	return
		h1.dna == h2.dna &&
		h1.start == h2.start &&
		h1.end == h2.end &&
		h1.strand == h2.strand;
}


// bool less_hit_index_ptr::operator()(HIT_P h1, HIT_P h2) const { return h1->index < h2->index; }

//the following semantics
//we want NULL < h and h !< NULL
//otherwise, they are ordered by score. if scores
//are equal, they are ordered by memory
//we want just a filler that allows us to not check the size of the hits, but
//


//
bool less_hit_score_ptr::operator()(HIT_PC h1, HIT_PC h2) const { 
	return 
		h1->score < h2->score ||
		(h1->score == h2->score && 
		 (h1->id < h2->id));
// 		 (h1->scoring_alg->id < h2->scoring_alg->id ||
// 		  (h1->scoring_alg->id == h2->scoring_alg->id &&
// 		   (h1->subname < h2->subname ||
// 			(h1->subname == h2->subname &&
// 			 (*(REG_P)h1 < *(REG_P)h2))))));

}






ostream& operator<<(ostream& f,HIT_CR h) {
	f<<FromDNAString(h.seq())<<'\t'
	 <<h.score<<'\t'
	 <<h.start<<'\t'
	 <<h.end<<'\t'
// 	 <<h.dna.species<<'\t'
// 	 <<h.dna.name<<'\t'
	 <<cis::PrintDNAStrand(h.strand);
		//<<h.scoring_alg->id;
// 	 <<PrintEnum(h.base.synteny_type())<<'\t'
// 	 <<PrintIsoforms(h.isoforms())<<'\t'
// 	 <<PrintGenes(h.genes);

 	return f;
}


  void PrintHitClusters(REGIONS_MULTI const& hits,
                        FILE * outstream){

    REGIONS_MULTI::iterator hits_iter;
    std::set<int>::iterator tags_iter;

    for (hits_iter = hits.begin(); hits_iter != hits.end(); ++hits_iter){
      HIT_CR h = *(HIT_P)(*hits_iter);
      for (tags_iter = h.cluster_tags.begin();
           tags_iter != h.cluster_tags.end(); ++tags_iter)
        fprintf(outstream, "%i\t%i\n", h.id, (*tags_iter));
    }

  }

  //append
  void PrintHits(REGIONS_MULTI const& hits, dna_t const& dna,
                 FILE * outstream){

    REGIONS_MULTI::iterator hits_iter;
    for (hits_iter = hits.begin(); hits_iter != hits.end(); ++hits_iter){
 
     HIT_CR h = *(HIT_P)(*hits_iter);
      cis::dnastring hitseq = h.strand == NEG ? RevCompDNA(h.seq()) : h.seq();

      fprintf(outstream, "%i\t%s\t%s\t%i\t%"PRId64"\t%"PRId64"\t%s\t%s\t%s\n",
              h.id, h.name().c_str(),
              FromDNAString(hitseq).c_str(), h.score, 
              h.start, h.end, dna.name.c_str(), dna.species().c_str(),
              cis::PrintDNAStrand(h.strand).c_str());
      
    }

  }

	

} // namespace cis
