#ifndef _MODULE_H
#define _MODULE_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "nuctrie.h"
#include "hit.h"
#include "permutation_collection.h"

namespace cis { 

  class hit; 
  class pattern; 

  
  class Module : public PermutationCollection {
    
  public:

    typedef std::map<std::string, hit_trie *> PATTERN_MAP;
    typedef PATTERN_MAP::const_iterator PATTERN_MAP_ITER;


    int id;
    cis::hit const** sites;
    cis_enum::rel_position relation_to_gene;
    int distance_to_gene;
    int _size;
    inline int size() const { return _size; };
    int64_t start;
    int64_t end;
    
  Module(int _id, cis::hit const** source_sites, int size) 
    : id(_id), _size(size) {
      sites = new cis::hit const*[size];
      for (int i = 0; i != size; ++i) sites[i] = source_sites[i];
      start = sites[0]->start;
      end = sites[_size - 1]->end;
    }
    
    ~Module(){
      delete sites; 
      //only deletes the array itself.  sites are owned by other...
    }
    

    Module MakeReverseComplement() const;
    
    //use this if the Module does indeed own the sites.
    void DeleteSites();

    //Bayesian probability that the module is transcriptionally functional
    float logFunctionality() const;
    
    static PATTERN_MAP patterns;

    static inline bool ValidElementPair(cis::hit const* query,
                                        cis::hit const* target){
      return query->name() == target->name();
    }


    static void Initialize(PATTERN_MAP const& patterns,
                           POINTS offset_curve,
                           POINTS permutation_curve,
                           POINTS orientation_curve,
                           POINTS hits_similarity_curve);

   
    

  };



} // end namespace cis


#endif // _MODULE_H
