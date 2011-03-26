#include "module.h"
#include "mapping.h"

#include <cmath>

#include "hit.h"
//#include "dna_types.h"

namespace cis {


  //make a copy of each source hit, in an order and orientation that 
  //represents the entire module reverse complemented
  //dest_module_hits must contain enough pointers, but this function
  //allocate the hits themselves
  Module Module::MakeReverseComplement() const {

    Module revcomp_module(this->id, this->sites, this->size());

    //determine the bounds of the source
    int64_t high_bound = this->sites[this->size() - 1]->end;

    for (int ind = 0; ind != this->size(); ++ind){
      cis::hit const & s = *this->sites[this->size() - 1 - ind];

      revcomp_module.sites[ind] = new 
        cis::hit(s.dna, cis::RevCompDNA(s.seq()), high_bound - s.end,
                 high_bound - s.start, ComplementStrand(s.strand),
                 s.score, s.id, s.name(), s.cluster_group);
    }

    return revcomp_module;
  }


  void Module::DeleteSites(){
    for (int ind = 0; ind != this->size(); ++ind) 
      delete this->sites[ind];
  }



  //calculate the log of the Bayesian probability P(functional|sequence)
  //the probability that the module is functional given its sequence
  float Module::logFunctionality() const {
    
    float log_probability_functional = 0.0;
    for (int i = 0; i != this->size(); ++i)
      log_probability_functional += this->sites[i]->score;

    return log_probability_functional;
  }


  void Module::Initialize
  (std::map<std::string, hit_trie *> const& _patterns,
   POINTS _offset_curve,
   POINTS _permutation_curve,
   POINTS _orientation_curve,
   POINTS _hits_similarity_curve){
    
    PermutationCollection::Initialize(_offset_curve,
                                      _permutation_curve,
                                      _orientation_curve,
                                      _hits_similarity_curve);
      
    patterns = _patterns;
    
  }

  
 

} //end namespace cis

