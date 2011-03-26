#include "regulon.h"
#include "mapping.h"

#include <cmath>

#include "module.h"


namespace cis {

  



  void Regulon::Initialize(std::vector<std::pair<cis_enum::rel_position,
                           cis_enum::rel_position> > const& allowed_pairings,
                           PermutationCollection::POINTS _offset_curve,
                           PermutationCollection::POINTS _permutation_curve,
                           PermutationCollection::POINTS _orientation_curve,
                           PermutationCollection::POINTS _hits_similarity_curve){
    
    PermutationCollection::Initialize(_offset_curve,
                                      _permutation_curve,
                                      _orientation_curve,
                                      _hits_similarity_curve);
    
    for (int i = 0; i < cis_enum::kNumRelPositions; ++i)
      for (int j = 0; j < cis_enum::kNumRelPositions; ++j)
        allowed_spatial_pairings[i][j] = false;
    
    for (size_t i = 0; i < allowed_pairings.size(); ++i)
      allowed_spatial_pairings
        [allowed_pairings[i].first]
        [allowed_pairings[i].second] = true;

  }

  
 
  
  bool Regulon::ValidElementPair(cis::Module const* query,
                                 cis::Module const* target){
    
    return Regulon::allowed_spatial_pairings
      [query->relation_to_gene]
      [target->relation_to_gene];

  }

  
  
  
} // end namespace cis


