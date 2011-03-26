#include <map>
#include <utility>
#include <string>


struct gene_relation {
  int relation_id;
  int cluster_id;
  cis_enum::rel_position relation;
  int offset;
  int gene_id;
  gene_relation(int _cluster_id, cis_enum::rel_position _relation,
                      int _offset, int _gene_id) :
    cluster_id(_cluster_id), relation(_relation),
    offset(_offset), gene_id(_gene_id){ }

};

typedef std::multimap<int, int> CLUSTER_MAP; //keyed on cluster_id, value = hit_id
typedef std::multimap<int, gene_relation *> RELATION_MAP; //keyed on gene_id
typedef std::map<int, cis::Module *> MODULE_MAP; //keyed on gene_id
typedef std::map<int, cis::Regulon *> REGULON_MAP; //keyed on gene_id
typedef std::map<int, cis::hit *> HITS_MAP; // key = hit_id
typedef std::multimap<int, int> ORTHOMAP; //ortholog_id, gene_id


MODULE_MAP BuildModules(CLUSTER_MAP const& clusters,
                        RELATION_MAP const& gene_relations);

cis::Regulon MakeRegulon(MODULE_MAP const& modules, int gene_id);


