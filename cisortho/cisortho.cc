//Program for assigning conserved functionality scores to genes
//based on prevalence of clustered binding sites

#include <getopt.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "enum.h"
#include "regulon.h"
#include "module.h"
#include "mapping.h"
#include "dna.h"
#include "conservation.h"

#include "cisortho_funcs.h"

/*
  input:

  hit_table                : hit_id dna start end strand seq score
  cluster_table            : cluster_id hit_id
  assoc_table              : relation_id cluster_id relation offset gene_id
  ortholog_table           : ortholog_id species gene_id

  output:
  conservation_table       : conservation_id ortholog_id gene_id support_gene_id conserved_reg_score
  module_function_table    : cluster_id log_prob_functionality
  mod_mapping_table        : conservation_id query_cluster_id target_cluster_id same_orientation log_prob_orthology
  hit_mapping_table        : query_cluster_id target_cluster_id query_hit_id target_hit_id log_prob_orthology


  Note, we do not even need gene exon coordinates for this.  All will be
  described by the relation and offset.

  And, what will user-friendly output look like?
  Do we also want a one-to-one site mapping within modules?
  That would be useful...

  In fact, there is really no reason to perform the selection of the top-scoring
  supporting gene in the actual c++ code.  Since all of the heavy lifting is already
  done, simply reporting the extra scores will not add to the major runtime.

*/




int main(int argc, char ** argv){

  char * clusters_file;
  char * gene_relation_file;
  char * hits_file;
  char * orthologs_file;
  char * guide_species;
  char * module_curve_file;
  char * regulon_curve_file;
  char * allowed_pairings_file;

  //output
  char * conservation_file;
  char * module_mapping_file;
  char * module_functionality_file;
  

  extern char *optarg;
  extern int optopt;
  char c;

  while ((c = getopt(argc, argv, "c:a:h:o:g:m:r:p:")) != -1) {
    switch(c) {
    case 'c': clusters_file = optarg; break;
    case 'a': gene_relation_file = optarg; break;
    case 'h': hits_file = optarg; break;
    case 'o': orthologs_file = optarg; break;
    case 'g': guide_species = optarg; break;
    case 'm': module_curve_file = optarg; break;
    case 'r': regulon_curve_file = optarg; break;
    case 'P': allowed_pairings_file = optarg; break;
    case '?': printf("Unknown argument %c\n", optopt); break;
    }
  }
  
  //parse allowed_spatial_pairings
  std::vector<std::pair<cis_enum::rel_position,
    cis_enum::rel_position> > allowed_spatial_pairings;

  

  


  //parse motif_patterns
  cis::Module::PATTERN_MAP motif_patterns;

  std::set<int> ortholog_ids;

  ORTHOMAP query_orthologs; //guide species orthologs
  ORTHOMAP target_orthologs; // non-guide species orthologs

  RELATION_MAP gene_relations;

  CLUSTER_MAP clusters;

  HITS_MAP hits;


  int hit_id;
  char hit_name[256];
  char hit_variant[256];
  char hit_sequence[4096];
  int hit_score;
  int64_t hit_start;
  int64_t hit_end;
  char hit_dna_name[256];
  char hit_species[256];
  char hit_strand[256];

  typedef std::set<cis::dna_t const*, cis::less_dna_ptr> DNAS;
  DNAS chromosomes;
  DNAS::const_iterator chrom_iter;
  std::pair<DNAS::iterator, bool> chrom_insert;


  //Load hits
  FILE * hits_file_p = fopen(hits_file, "r");
  while (! feof(hits_file_p)){
    fscanf(hits_file_p, "%i\t%s\t%s\t%s\t%i\t%"PRId64"\t%"PRId64"\t%s\t%s\t%s\n",
           &hit_id, hit_name, hit_variant, hit_sequence, &hit_score,
           &hit_start, &hit_end, hit_dna_name, hit_species, hit_strand);
    
    cis::dna_t const dna("", hit_species, hit_dna_name);
    chrom_iter = chromosomes.find(&dna);
    if (chrom_iter == chromosomes.end())
      chrom_insert = 
        chromosomes.insert(new cis::dna_t("", hit_species, hit_dna_name));

    cis::dna_t const& this_dna = *(*chrom_insert.first);
    hits[hit_id] = 
      new cis::hit(this_dna, 
                   cis::dnastring(std::string(hit_sequence)),
                   hit_start, hit_end, 
                   (cis::dna_strand)cis_enum::from_string(std::string(hit_strand)),
                   hit_score, hit_id, hit_name);
  }


  fclose(hits_file_p);

  //Load clusters
  int cluster_id;
  FILE * clusters_file_p = fopen(clusters_file, "r");
  while (! feof(clusters_file_p)){
    fscanf(clusters_file_p, "%i\t%i\n", &cluster_id, &hit_id);
    clusters.insert(std::pair<int,int>(cluster_id, hit_id));
  }
  fclose(clusters_file_p);


  //Load gene_relations
  int cluster_offset;
  char gene_rel_position[256];
  int gene_id;

  FILE * gene_relation_file_p = fopen(gene_relation_file, "r");
  while (! feof(gene_relation_file_p)){
    fscanf(gene_relation_file_p, "%i\t%s\t%i\t%i\n", 
           &cluster_id, gene_rel_position, &cluster_offset,
           &gene_id);
    
    cis_enum::rel_position relation =
      static_cast<cis_enum::rel_position>
      (cis_enum::from_string(std::string(gene_rel_position)));

    gene_relations.insert
      (std::pair<int, gene_relation *>
       (gene_id, new gene_relation(cluster_id, relation, 
                                         cluster_offset, gene_id)));
    
  }
  fclose(gene_relation_file_p);


  //Load orthologs and ortholog_ids
  int ortholog_id;
  char species[256];

  FILE * orthologs_file_p = fopen(orthologs_file, "r");
  while (! feof(orthologs_file_p)){
    fscanf(orthologs_file_p, "%i\t%s\t%i\n",
           &ortholog_id, species, &gene_id);

    ortholog_ids.insert(ortholog_id);

    if (strcmp(species, guide_species) == 0)
      query_orthologs.insert(std::pair<int, int>(ortholog_id, gene_id));
    else
      target_orthologs.insert(std::pair<int, int>(ortholog_id, gene_id));
                            
  }
  fclose(orthologs_file_p);


  
    

  //build all modules.  somewhat wasteful, but oh well...
  //modules correspond with clusters, will be keyed by cluster_id
  MODULE_MAP modules = BuildModules(clusters, gene_relations);

  FILE * conservation_file_p = fopen(conservation_file, "w");
  FILE * module_mapping_file_p = fopen(module_mapping_file, "w");
  FILE * module_functionality_file_p = fopen(module_functionality_file, "w");

  int conservation_id = 0;


  std::vector<ComponentScoreData> module_counts;
  IterateOrthologPairs(ortholog_ids, query_orthologs, target_orthologs,
                       regulons, &EstimateModuleCounts, &module_counts);

  //build module stats curves, initialize Module::
  //parse module curve file
  cis::PermutationCollection::POINTS
    module_offset_curve,
    module_permutation_curve,
    module_orientation_curve,
    module_hits_similarity_curve;

  cis::Module::Initialize(motif_patterns,
                          module_offset_curve,
                          module_permutation_curve,
                          module_orientation_curve,
                          module_hits_similarity_curve);


  std::vector<ComponentScoreData> regulon_counts;
  IterateOrthologPairs(ortholog_ids, query_orthologs, target_orthologs,
                       regulons, &EstimateRegulonCounts, &regulon_counts);

  //parse regulon curve file
  cis::PermutationCollection::POINTS
    regulon_offset_curve,
    regulon_permutation_curve,
    regulon_orientation_curve,
    regulon_hits_similarity_curve;

    
  //Initialize Regulon and Module classes
  cis::Regulon::Initialize(allowed_spatial_pairings,
                           regulon_offset_curve,
                           regulon_permutation_curve,
                           regulon_orientation_curve,
                           regulon_hits_similarity_curve);


  //iterate through ortholog_group, query ortholog...
  //find the best N supports for each species, and print everything out
  //in a special format
  std::set<int>::iterator ortholog_iter;
  ORTHOMAP::iterator query_ortho_iter, target_ortho_iter;
  std::pair<ORTHOMAP::iterator, ORTHOMAP::iterator> query_range, target_range;

  for (ortholog_iter = ortholog_ids.begin();
       ortholog_iter != ortholog_ids.end(); ++ortholog_iter){
    int ortholog_id = *ortholog_iter;
    
    query_range = query_orthologs.equal_range(ortholog_id);
    
    for (query_ortho_iter = query_range.first;
         query_ortho_iter != query_range.second;
         ++query_ortho_iter){
      
      cis::Regulon const* query = regulons.find((*query_ortho_iter).second);
      
      OrthologSupport query_support =
        BestSupportingOrthologs(ortholog_id, max_supporting_orthologs,
                                query, target_orthologs, regulons);

      PrintConservedOrthologs(conservation_file_p, query_support);
      
    }
  }

    


  //cycles through each ortholog group, then
  //each set of query orthologs (guide species)
  //each set of target orthologs, constructs the regulon from
  //the set of modules (indexed by gene_id)

    
    /*
    //print conservation table
    fprintf(conservation_file_p, "%i\t%i\t%i\t%i\t%f\n",
            conservation_id, ortholog_id, query_gene_id,
                target_gene_id, conserved_reg_score);
      

        for (int mm = 0; mm != query_regulon.size(); ++mm){

          int target_ind = module_mapping[mm].best_target_index;
          if (target_ind != -1){

            //print module mapping table
            fprintf(module_mapping_file_p, "%i\t%i\t%i\t%s\t%f\n",
                    conservation_id, query_regulon.sites[mm]->id,
                    target_regulon.sites[target_ind]->id,
                    (module_mapping[mm].same_orientation ? "same" : "opposite"),
                    module_mapping[mm].log_orthology);
            
            //print module functionality
            fprintf(module_functionality_file_p, "%i\t%f\n",
                    query_regulon.sites[mm]->id,
                    query_regulon.sites[mm]->logFunctionality());

            fprintf(module_functionality_file_p, "%i\t%f\n",
                    target_regulon.sites[target_ind]->id,
                    target_regulon.sites[target_ind]->logFunctionality());
           

          } // if
        } // modules

        delete [] module_mapping;

      } // target ortholog
    } // query ortholog
  } // ortholog group
    */

  fclose(conservation_file_p);
  fclose(module_mapping_file_p);
  fclose(module_functionality_file_p);

  //clean up

} // main
