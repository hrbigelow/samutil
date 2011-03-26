//Functions for scoring conservation of hits across genomes
#include "conservation.h"
#include "mapping.h"
#include "regulon.h"
#include "module.h"

#include <cmath>

namespace cis {

  //computes regulon similarity, storing the module mapping pattern
  float RegulonOrthologyScore(cis::Regulon const& query,
                              cis::Regulon const& target,
                              MappingData * best_module_mapping){
    
    float summary_score = -HUGE_VAL;
    MappingData * module_mapping = new MappingData[query.size()];
    extend_mapping(query, target, 0, module_mapping, best_module_mapping,
                   &summary_score);
    delete module_mapping;
    return summary_score;
  }


  float ModuleOrthologyScore(cis::Module const& query,
                             cis::Module const& target,
                             MappingData * best_mapping,
                             bool * same_orientation){

    
    float forward_score = -HUGE_VAL;

    MappingData * forward_mapping = new MappingData[query.size()];
    extend_mapping(query, target, 0, forward_mapping, 
                   best_mapping, &forward_score);
    
    float reverse_score = forward_score;

    MappingData * reverse_mapping = new MappingData[query.size()];
    cis::Module target_revcomp = target.MakeReverseComplement();
    extend_mapping(query, target_revcomp, 0, reverse_mapping, 
                   best_mapping, &reverse_score);

    *same_orientation = forward_score > reverse_score;

    MappingData * selected_mapping = *same_orientation ?
      forward_mapping : reverse_mapping;

    //copy selected_mapping to site_mapping
    for (int qi = 0; qi != query.size(); ++qi)
      best_mapping[qi] = selected_mapping[qi];

    delete forward_mapping;
    delete reverse_mapping;
    
    return *same_orientation ? forward_score : reverse_score;

  }


  //process the query and target, saving the best N top genes in each supporting
  //species for each given query.  
  class SupportGene {

    typedef std::vector<std::pair<MappingData *, float> > MODULE_PAIRS;
    
    cis::Regulon & target;
    float functional_orthology_score;
    MappingData * regulon_mapping;
    MODULE_PAIRS module_mappings;
    

    ~SupportGene(){
      delete regulon_mapping;
      for (size_t mi = 0; mi != module_mappings.size(); ++mi)
        if (module_mappings[mi].first != NULL) {
          delete module_mappings[mi].first;
          module_mappings[mi].first = NULL;
        }
    }
  };
  
  
  struct less_support_gene {
    bool operator()(SupportGene const& a,
                    SupportGene const& b) const {
      return
        a.functional_orthology_score <
        b.functional_orthology_score;
    }
  };


  class OrthologSupport {
    
    cis::Regulon * query;
    typedef std::set<SupportGene, less_support_gene> SUPPORT_GENES;
    typedef std::map<std::string, SUPPORT_GENES> SUPPORT_SPECIES;
    
    SUPPORT_SPECIES support_species;
    
  };


  //
  std::string DNADifference(cis::dnastring const& query_seq,
                            cis::dnastring const& target_seq,
                            bool same_orientation){

    if (query_seq.size() != target_seq.size()){
      printf(stderr, "DNADifference: different sized query and target sequences.");
      exit(50);
    }

    std::string difference(target_seq.size());

    if (same_orientation){
      for (size_t qi = 0; qi != query_seq.size(); ++qi)
        difference[qi] = query_seq[qi] == target_seq[qi] ? 
          '.' : Nuc2char(target_seq[qi]);
    } else {
      for (size_t qi = 0; qi != query_seq.size(); ++qi)
        difference[qi] = query_seq[qi] == target_seq[qi] ? 
          ',' : Nuc2lcchar(target_seq[qi]);
    }
    return difference;

  }

 

  //print the pattern of conserved hits within modules, using
  //a symbolic string...
  //for each pair of mapped modules, print the entire
  //combinatorial alignment below the original query
  //module...
  void PrintConservedOrthologs(FILE * conservation_file,
                               OrthologSupport conserved_gene){
    for (int qi = 0; qi != conserved_gene.query->size(); ++qi){
      cis::Module const& qmodule = *conserved_gene.query[qi];
      fprintf(conservation_file, "%i\t%3.5f\t%i\t%s\t",
              qmodule.id, 10000.0,
              qmodule.distance_to_gene,
              cis_enum::print(qmodule.relation_to_gene).c_str());

      for (int hi = 0; hi != qmodule.size(); ++hi){
        if (hi != 0) {
          fprintf(conservation_file, "---|%-5i|---", 
                  qmodule.sites[hi]->start -
                  qmodule.sites[hi-1]->end);
        }
        fprintf(conservation_file, "%s", 
                FromDNAString(qmodule.sites[hi].seq()).c_str());
      }
      //now print every matched module...
      OrthologSupport::SUPPORT_SPECIES & sspec = conserved_gene.support_species;
      OrthologSupport::SUPPORT_SPECIES::iterator s_iter;

      for (s_iter = sspec.begin(); s_iter != sspec.end(); ++s_iter){
        std::string & target_species = (*s_iter).first;
        OrthologSupport::SUPPORT_GENES & sgenes = (*s_iter).second;
        OrthologSupport::SUPPORT_GENES::iterator sgenes_iter;
        for (sgenes_iter = sgenes.begin(); 
             sgenes_iter != sgenes.end(); ++sgenes_iter){

          SupportGene & sgene = (*sgenes_iter);
          fprintf(conservation_file, "%i\t%3.5f\t%i\t%s\t",
                  sgene.target.id, sgene.functional_orthology_score,
                  sgene.target.distance_to_gene,
                  cis_enum::print(sgene.target.relation_to_gene).c_str());
          
          for (int hi = 0; hi != sgene.target.size(); ++hi){
            if (hi != 0) {
              fprintf(conservation_file, "---|%-5i|---", 
                      
                      sgene.target.sites[hi]->start -
                      sgene.target.sites[hi-1]->end);
            }
            fprintf(conservation_file, "%s", 
                    FromDNAString(sgene.target.sites[hi].seq()).c_str());
          }
          
        }


    
  }

  
  //call this for each query ortholog.  It calculates at most the top N
  //supporting orthologs in each non-guide species, and prints them out
  //in a summary format.
  OrthologSupport BestSupportingOrthologs(int ortholog_id, 
                                          int max_supporting_orthologs,
                                          cis::Regulon const* query,
                                          ORTHOMAP const& target_orthologs,
                                          REGULON_MAP const& regulons){
    
    ORTHOMAP::iterator target_ortho_iter;
    std::pair<ORTHOMAP::iterator, ORTHOMAP::iterator> target_range;
    
    OrthologSupport ortho_support(max_supporting_orthologs);
    
    target_range = target_orthologs.equal_range(ortholog_id);

    for (target_ortho_iter = target_range.first;
         target_ortho_iter != target_range.second;
         ++target_ortho_iter){
    
      cis::Regulon const* target = regulons.find((*target_ortho_iter).second);
      
      //try out this target against the query
      UpdateQuerySupport(*query, *target, target_species,
                         max_supporting_orthologs,
                         &ortholog_support);
    }

    //now calculate the module submappings for each qualifying ortholog
    //assume the support genes sets are the proper size,
    //and create module submappings for each
    OrthologSupport::SUPPORT_SPECIES::iterator species_iter;
    OrthologSupport::SUPPORT_GENES::iterator genes_iter;
    for (species_iter = ortho_support.support_species.begin();
         species_iter != ortho_support.support_species.end();
         ++species_iter){
      for (genes_iter = (*species_iter).second.begin();
           genes_iter != (*species_iter).second.end();
           ++genes_iter){
        SupportGene & sgene = (*genes_iter);
        sgene.module_mappings = 
          GetModuleSubmappings(*query, *sgene.target, sgene.regulon_mapping);
      }      
      
    }
    
    return ortho_support;
  }


  //updates the supporting orthologs.
  //calls 
  void UpdateQuerySupport(cis::Regulon const& query,
                          cis::Regulon const& target,
                          std::string const& target_species,
                          int max_ortholog_support,
                          OrthologSupport * ortholog_support){
    
    MappingData * best_module_mapping = new MappingData[query.size()];
    float orthology_score = 
      RegulonOrthologyScore(query, target, best_module_mapping);

    OrthologSupport::SUPPORT_SPECIES::iterator
      support_iter = ortholog_support->support_species.find(target_species);
    
    if (support_iter == ortholog_support->end()){
      std::pair<OrthologSupport::SUPPORT_SPECIES::iterator, bool> support_insert =
        ortholog_support->insert(std::make_pair(target_species, 
                                                OrthologSupport::SUPPORT_GENES()));
      support_iter = support_insert.first;
    }

    OrthologSupport::SUPPORT_GENES & sgenes = (*support_iter).second;
    
    sgenes.insert(SupportGene(target, orthology_score, best_module_mapping, 
                              SupportGene::MODULE_PAIRS));
    
    while (sgenes.size() > max_ortholog_support)
      sgenes.erase(sgenes.rbegin());
    

  }
    

  

                                  
  //returns the element wise mappings for each module / module pair
  //in module_mapping
      SupportGene::MODULE_PAIRS
      GetModuleSubMappings(cis::Regulon const& query,
                           cis::Regulon const& target,
                       MappingData const* module_mapping){

    std::vector<MappingData *> submappings(query.size());
    std::vector<std::pair<MappingData *, float> > 
      best_submappings(query.size());

    for (int qi = 0; qi != query.size(); ++qi){
      submappings[qi] = new MappingData[query.sites[qi]->size()];
      best_submappings[qi] = 
        std::make_pair(new MappingData[query.sites[qi]->size()], -HUGE_VAL);

      if (module_mapping[qi].target_used)
        extend_mapping(*query.sites[qi],
                       *target.sites[module_mapping[qi].target_index],
                       0, submappings[qi], best_submappings[qi].first,
                       &best_submappings[qi].second);
    }
    for (int qi = 0; qi != query.size(); ++qi)
      delete submappings[qi];

    return best_submappings;

  }

    
  //For accumulating best component Module statistics.

  //Here we assume no CPD curves are available for Module or Regulon.
  //So call extend_mapping on the query and target module, computing
  //the best component-wise mapping
  float OrthologyScore(cis::Module const* query,
                       cis::Module const* target,
                       bool * same_orientation){
    
    MappingData * mapping = new MappingData[query->size()];
    
    float score = ModuleOrthologyScore(*query, *target, mapping, 
                                       same_orientation);

    delete mapping;
    return score;
  }  



  //calculates the average fluctuation in the score in the sequence of site
  //variants from the query to the target, one point mutation at a time.
  float OrthologyScore(cis::hit const* query_hit, 
                       cis::hit const* target_hit,
                       bool * same_orientation){
    
    assert(query_hit->name() == target_hit->name()); 
    //fails if we have terminal_id, must fix

    assert(query_hit->seq().size() == target_hit->seq().size());

    Module::PATTERN_MAP_ITER pattern_iter = 
      Module::patterns.find(query_hit->name());

    if (pattern_iter == Module::patterns.end()){
      fprintf(stderr, 
              "MotifAvgScoreChange: Couldn't find the hit pattern called %s.  Returning -1",
              query_hit->name().c_str());

      return -1.0;
    } else {
      //now generate the series scores, starting from the query site and ending
      //with the target site, each time, advancing to the next position where
      //they differ.
      hit_trie * patnode = (*pattern_iter).second;
      cis::dnastring current_sequence = query_hit->seq();
      cis::dnastring target_sequence = target_hit->seq();

      int site_length = query_hit->seq().size();
      int * score_path = new int[site_length];
      score_path[0] = query_hit->score;

      //populate score_path
      for (int current_index = 1; current_index != site_length; ++current_index){
        if (current_sequence[current_index] == target_sequence[current_index]){
          //no change, just use last score
          score_path[current_index] = score_path[current_index-1];
        } else {
          current_sequence[current_index] = target_sequence[current_index];
          score_path[current_index] = patnode->score(current_sequence);
        }
      }
    
      //calculate average absolute value of score change
      float sum_of_changes = 0;

      for (int i = 1; i != site_length; ++i)
        sum_of_changes += abs(score_path[i] - score_path[i-1]);

      float average_score_change = 
        sum_of_changes / static_cast<float>(site_length);

      delete score_path;
    
      return average_score_change;
    
    } 

  }



  //go through every legal pairing of query and target modules,
  //and calculate the best component score counts
  void EstimateModuleCounts(cis::Regulon const& query,
                            cis::Regulon const& target,
                            std::vector<ComponentScoreData> * counts){

    int target_size = target.size();

    cis::Module ** target_revcomp = new cis::Module *[target_size];
    int target_last_ind = target_size - 1;
    for (int ti = 0; ti != target_size; ++ti)
      target_revcomp[ti] = 
        new cis::Module(target.sites[target_last_ind - ti]->MakeReverseComplement());

    ComponentScoreData best_scores(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL);

    int twice_target_size = target_size * 2;
    cis::Module const** all_target_sites = new cis::Module const*[target_size * 2];
    for (int ti = 0; ti != target_size; ++ti)
      all_target_sites[ti] = target.sites[ti];
    for (int ti = target_size; ti != twice_target_size; ++ti)
      all_target_sites[ti] = target_revcomp[ti - target_size];


    //compare each query with each target site both forward and
    //reverse complemented
    for (int qi = 0; qi != query.size(); ++qi){
      MappingData * mapping = new MappingData[query.sites[qi]->size()];
      MappingData * best_mapping = new MappingData[query.sites[qi]->size()];
      for (int ti = 0; ti != twice_target_size; ++ti){
        if (Regulon::ValidElementPair(query.sites[qi], all_target_sites[ti])){
          extend_mapping(*query.sites[qi], *all_target_sites[ti], 0,
                         mapping, best_mapping, &best_scores);
          counts->push_back(best_scores);
        }
        
        //zero the mapping in between each comparison
        for (int qmi = 0; qmi != query.sites[qi]->size(); ++qmi)
          mapping[qmi].target_used = false;
        
      }
      delete mapping;
      delete best_mapping;
    }

    //clean up
    for (int ti = 0; ti != target_size; ++ti) delete target_revcomp[ti];
    delete target_revcomp;
    delete all_target_sites;

  }


  //augment counts with one tetrad of best score components
  void EstimateRegulonCounts(cis::Regulon const& query,
                             cis::Regulon const& target,
                             std::vector<ComponentScoreData> * counts){

    ComponentScoreData best_scores(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
    MappingData * mapping = new MappingData[query.size()];
    MappingData * best_mapping = new MappingData[query.size()];

    extend_mapping(query, target, 0, mapping, best_mapping, &best_scores);
    delete mapping;
    delete best_mapping;
    counts->push_back(best_scores);

  }


  

  //return the number of mapped elements in mapping
  int NumberMapped(MappingData * mapping, int query_size){

    int number_mapped = 0;
    for (int qi = 0; qi < query_size; ++qi)
      if (mapping[qi].target_index != -1) 
        number_mapped++;

    return number_mapped;
  }


  //computes the current component-wise scores from the mapping data
  ComponentScoreData ScoreSummary(MappingData * mapping, int query_size){
    
    int number_mapped = NumberMapped(mapping, query_size);

    //first, calculate the average position
    int64_t total_offset = 0;

    ComponentScoreData score(0.0,0.0,0.0,0.0);
    
    for (int qi = 0; qi < query_size; ++qi)
      if (mapping[qi].target_used)
        total_offset += mapping[qi].offset;
    
    int64_t average_offset = total_offset / number_mapped;
    
    //adjust query_target_offset
    for (int qi = 0; qi < query_size; ++qi)
      if (mapping[qi].target_used){
        score.offset_score += 
          abs(mapping[qi].offset - average_offset);

        if (mapping[qi].in_order) score.permutation_score++;
        if (mapping[qi].same_orientation) score.orientation_score++;
        score.similarity_score += mapping[qi].similarity;
      }

    return score;

  }


  //assume the mapping 




  //use log-sum-exp trick to get a sum
  //of potentially undeflow numbers of all mapped modules.
  /*
  float ConservedRegulationScore(cis::Regulon const& query,
                                 cis::Regulon const& target,
                                 MappingData * mapping){

    float min_log_functional_orthology = HUGE_VAL;

    float *log_functional_orthology = new float[query.size()];
    float *scaled_raw_probability = new float[query.size()];
    
    for (int qi = 0; qi < query.size(); ++qi){
      if (mapping[qi].best_target_index != -1){
        log_functional_orthology[qi] =
          query.sites[qi]->logFunctionality() +
          target.sites[mapping[qi].best_target_index]->logFunctionality() +
          ModuleOrthologyScore(*query.sites[qi], 
                               *target.sites[mapping[qi].best_target_index],
                               site_mapping,
                               same_orientation);

        if (mapping[qi].log_orthology_score < min_log_functional_orthology)
          min_log_functional_orthology = mapping[qi].log_orthology_score;
      }
    }

    for (int qi = 0; qi < query.size(); ++qi){
      if (mapping[qi].best_target_index != -1){
        scaled_raw_probability[qi] = 
          exp(log_functional_orthology[qi] - min_log_functional_orthology);
      }
    }

    float sum_scaled_raw_prob = 0.0;
    for (int qi = 0; qi < query.size(); ++qi){
      if (mapping[qi].best_target_index != -1){
        sum_scaled_raw_prob += scaled_raw_probability[qi];
      }
    }

    float log_sum_functional_orthology =
      log(sum_scaled_raw_prob) + min_log_functional_orthology;

    return log_sum_functional_orthology;
  }
  */


    
} // end namespace cis



//Input:  Ortholog group and guide gene within it
//Output: A sorted list of one gene per non-guide species (called a support gene).
//Each selected support gene is the gene in the other species having the best
//GeneRegulonSimilarity among all of the genes in the support species.
//The genes are sorted descending by their GeneRegulonSimilarity and may be used
//as a 2-d ranking scheme (see below)


/*
  Best-rank.  For the entire set of genes in a given guide genome,
  assume you have for each gene a set of N pairwise
  GeneRegulonSimilarity scores between that gene and supporting genes.
  For S supporting species overall, one can produce S sorted lists of
  genes based on the GRS score of the s-th supporting gene.  For each
  guide gene then, the best rank obtained in any of these lists is its
  best-rank score.

  For genes lacking one or more of the S maximum possible supporting
  genes, the convention is that the GRS score is -inf.  This way, genes
  with no support are at the bottom of the list
*/
