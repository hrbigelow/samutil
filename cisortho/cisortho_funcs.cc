#include "cis_phylo_funcs.h"


//build a set of modules, keyed on the relation_id
MODULE_MAP BuildModules(CLUSTER_MAP const& clusters,
                        RELATION_MAP const& gene_relations){

  MODULE_MAP::iterator module_iter;
  std::pair<CLUSTER_MAP::iterator, CLUSTER_MAP::iterator> cluster_range;
  RELATION_MAP::iterator relation_iter;

  for (relation_iter = gene_relations.begin();
       relation_iter != gene_relations.end(); ++relation_iter){
    
    int relation_id = (*relation_iter).first;

    cluster_range = clusters.equal_range(cluster_id);
    int size = std::distance(cluster_range.first, cluster_range.second);
    cis::hit const ** sites = new cis::hit const * [size];
    int ci = 0;
    while (cluster_range.first != cluster_range.second)
      sites[ci++] = hits[(*cluster_range.first++)];
    
    modules[gene_id] = new cis::Module(relation_id, sites, size);

    delete sites;
    
  }
  
  return modules;
}




cis::Regulon MakeRegulon(MODULE_MAP const& modules, int gene_id){

  std::pair<MODULE_MAP::const_iterator, MODULE_MAP::const_iterator>
    module_range = modules.equal_range(gene_id);

  MODULE_MAP::const_iterator module_iter = module_range.first;

  int regulon_size = 
    std::distance(module_range.first, module_range.second);
  
  cis::Module const ** regulon_elems =
    new cis::Module const * [regulon_size];
  
  int mi = 0;
  while (module_iter != module_range.second)
    regulon_elems[mi++] = (*module_iter++).second;

  cis::Regulon regulon(regulon_elems, regulon_size);

  delete regulon_elems;

  return regulon;
  
}


template <typename Accumumator>
void IterateOrthologPairs(std::set<int> const& ortholog_ids,
                          ORTHOMAP const& query_orthologs,
                          ORTHOMAP const& target_orthologs,
                          REGULON_MAP const& regulons,
                          void (* Comparison)(cis::Regulon const& query,
                                              cis::Regulon const& target,
                                              Accumulator * accu),
                          Accumulator * accumulator){
  std::set<int>::iterator ortholog_iter;
  ORTHOMAP::iterator query_ortho_iter, target_ortho_iter;
  std::pair<ORTHOMAP::iterator, ORTHOMAP::iterator> query_range, target_range;

  for (ortholog_iter = ortholog_ids.begin();
       ortholog_iter != ortholog_ids.end(); ++ortholog_iter){
    int ortholog_id = *ortholog_iter;
    
    IterateOrthologPair(ortholog_id, query_orthologs, target_orthologs, 
                        regulons, Comparison, accumulator);
  }
}





template <typename Accumumator>
void IterateOrthologPair(int ortholog_id,
                         ORTHOMAP const& query_orthologs,
                         ORTHOMAP const& target_orthologs,
                         REGULON_MAP const& regulons,
                         void (* Comparison)(cis::Regulon const& query,
                                             cis::Regulon const& target,
                                             Accumulator * accu),
                         Accumulator * accumulator){
  ORTHOMAP::iterator query_ortho_iter, target_ortho_iter;
  std::pair<ORTHOMAP::iterator, ORTHOMAP::iterator> query_range, target_range;
    
  query_range = query_orthologs.equal_range(ortholog_id);
  target_range = target_orthologs.equal_range(ortholog_id);

  for (query_ortho_iter = query_range.first;
       query_ortho_iter != query_range.second;
       ++query_ortho_iter){

    cis::Regulon const* query = regulons.find((*query_ortho_iter).second);

    for (target_ortho_iter = target_range.first;
         target_ortho_iter != target_range.second;
         ++target_ortho_iter){
        
      cis::Regulon const* target = regulons.find((*target_ortho_iter).second);

      Comparison(*query, *target, accumulator);
    }
  }
}



//need a function for iterating through all guide/non-guide gene pairs,
//constructing a pair of regulons, and performing the proper
//pairwise comparison.
//there will be two uses of this.  The first is to build component score
//distributions, four each for module-level and regulon-level characteristics.

//After these statistics are smoothed into CPDs they are used to initialize
//the module and regulon classes, and a second round of comparisons are made,
//this time calculating an aggregate score using the CPDs for significance
//estimation.

/*TODO:  
  2.  distribution smoothing function
  3.  regulon and module member functions for ???
 */
