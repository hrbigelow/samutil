//Functions for scoring conservation of hits across genomes
#ifndef _CONSERVATION_H
#define _CONSERVATION_H

#include <vector>
#include <utility>

class MappingData;
class ComponentScoreData;

namespace cis {
  
  class Regulon;
  class Module;
  class hit;

  float RegulonOrthologyScore(cis::Regulon const& query,
                              cis::Regulon const& target,
                              MappingData * module_mapping);

  float ModuleOrthologyScore(cis::Module const& query,
                             cis::Module const& target,
                             MappingData * site_mapping,
                             bool * same_orientation);

  std::vector<std::pair<MappingData *, float> >
  GetModuleSubMappings(cis::Regulon const& query,
                       cis::Regulon const& target,
                       MappingData const* module_mapping);


  float OrthologyScore(cis::Module const* query,
                       cis::Module const* target,
                       bool * same_orientation);


  float OrthologyScore(cis::hit const* query_hit, 
                       cis::hit const* target_hit,
                       bool * same_orientation);

  void EstimateModuleCounts(cis::Regulon const& query,
                            cis::Regulon const& target,
                            std::vector<ComponentScoreData> * counts);

  void EstimateRegulonCounts(cis::Regulon const& query,
                             cis::Regulon const& target,
                             std::vector<ComponentScoreData> * counts);

  ComponentScoreData ScoreSummary(MappingData * mapping, int query_size);



  float ConservedRegulationScore(cis::Regulon const& query,
                                 cis::Regulon const& target,
                                 MappingData * mapping);
  
  
} // end namespace cis


//Input:  Ortholog group and guide gene within it
//Output: A sorted list of one gene per non-guide species (called a support gene).
//Each selected support gene is the gene in the other species having the best
//GeneRegulonSimilarity among all of the genes in the support species.
//The genes are sorted descending by their GeneRegulonSimilarity and may be used
//as a 2-d ranking scheme (see below)
//float WholeGeneScore(cis::ortholog_group const& orthologs);

#endif // _CONSERVATION_H
