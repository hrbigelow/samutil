#ifndef REGULON_H
#define REGULON_H

#include "permutation_collection.h"

#include <vector>
#include <map>
#include <string>

#include "enum.h"

class MappingData;

namespace cis {

  
  class Module;

  class Regulon : public PermutationCollection {

    static bool allowed_spatial_pairings
      [cis_enum::kNumRelPositions]
      [cis_enum::kNumRelPositions];
    
  public:

    cis::Module const** sites;
    int _size;
    inline int size() const { return _size; };
    
    //Regulons do not own their sites, only the pointers
  Regulon(cis::Module const** _sites, int size) :
    _size(size) { 
      sites = new cis::Module const *[_size];
      for (int i = 0; i != _size; ++i) sites[i] = _sites[i];
    }
    
    ~Regulon(){ delete sites; }
    
    
    //initialize a 2D vector
    static void Initialize(std::vector<std::pair<cis_enum::rel_position,
                           cis_enum::rel_position> > const& allowed_pairings,
                           POINTS offset_curve,
                           POINTS permutation_curve,
                           POINTS orientation_curve,
                           POINTS hits_similarity_curve);
      
    static bool ValidElementPair(cis::Module const* query,
                                 cis::Module const* target);


  };

} // end namespace cis

#endif //REGULON_H
  
