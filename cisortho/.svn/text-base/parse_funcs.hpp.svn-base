#ifndef _PARSE_FUNCS_HPP
#define _PARSE_FUNCS_HPP

#include <climits>
#include <vector>
#include <string>
#include <map>

#include "pwm.h"
#include "pattern.h"
#include "diversity_cluster.h"
#include "specific_cluster.h"
#include "dna_types.h"
// #include "sql.h"

#include "input_grammar.hpp"

using std::string;

namespace cis { 
//   struct sql_credentials; 
  class pattern;
}


//make a pattern
struct make_pattern {

  make_pattern(std::vector<cis::pattern *> & patterns_) : patterns(patterns_){}
  template <typename TreeT, typename IterT>
  void operator()(TreeT const& tree,
                  IterT const& /*first*/, IterT const& /*last*/) const;
  std::vector<cis::pattern *> & patterns;
};



struct make_cluster {
  make_cluster(std::vector<cis::cluster *> & clusters_) : clusters(clusters_){}
  template <typename TreeT, typename IterT>
  void operator()(TreeT const& tree,
                  IterT const& /*first*/, IterT const& /*last*/) const;
  std::vector<cis::cluster *> & clusters;
};


//set a bool from lowercase 'true' or 'false' string
struct assign_bool {
  
  assign_bool(bool & bool_val_) : bool_val(bool_val_){}
  template <typename IterT>
  void operator()(IterT const& first, IterT const& last) const;
  bool & bool_val;
};




struct parsed {

//   sql::sql_credentials sql_info;

  std::vector<cis::pattern *> patterns;
  make_pattern pattern_maker;
  
  std::vector<cis::cluster *> clusters;
  make_cluster cluster_maker;

  std::vector<std::string> species;

  assign_bool set_merge_overlap;
  bool merge_overlap;

  std::string data_directory;

  std::string posfile;
  std::string treefile;
  std::string dnaindexfile;
//   std::string sql_output_prefix;
  std::string hits_file;
  std::string clusters_file;

  int max_total_hits;
  int max_collection_size;
  int max_pwm_length;
  int min_pwm_length;
  int min_collection_score;
  int max_patterns;


  parsed() : pattern_maker(patterns),
             cluster_maker(clusters),
             set_merge_overlap(this->merge_overlap)
  {}

  void prepend_paths() {
    data_directory += "/";
    posfile = data_directory + posfile;
    treefile = data_directory + treefile;
    dnaindexfile = data_directory + dnaindexfile;
  }


};


//parse a vert_matrix, as defined by grammar rule 'vert_matrix'
//create an [N][4] matrix
template <typename TreeT>
IMATRIX vert_matrix(TreeT const& tree){
  
  IMATRIX result(12);
  std::vector<std::map<char, int> > result_map; //holds the intermediate matrix structure
  std::map<char, int>::iterator mit;

  //the iterator TreeT is the type of pointer to a tree node.
  //the siblings of a given node 't' are in the range [t->children.begin(), t->children.end())

  const size_t num_nucs = 4;

  //std::map<char, int> nucleotides;
  std::vector<char> nucs_order(num_nucs);

  TreeT header = tree.children[0];
  TreeT data = tree.children[1];
  TreeT data_line = data.children[0];
  size_t matrix_length = data.children.size();


  //here, rows are the sequence positions, columns are the four nucleotides
  for (size_t col = 0; col < num_nucs; ++col) 
    nucs_order[col] = *(header.children[col].value.begin());

  result_map.resize(matrix_length);
  result.resize(matrix_length);

  for (size_t row = 0; row < matrix_length; ++row){
    TreeT const& rownode = data.children[row];
    for (size_t col = 0; col < num_nucs; ++col){
      TreeT const& cell = rownode.children[col];
      std::string txt(cell.value.begin(), cell.value.end());
      
      result_map[row][nucs_order[col]] = atoi(txt.c_str());
    }
  }

  for (size_t row = 0; row < matrix_length; ++row){
    for (mit = result_map[row].begin(); mit != result_map[row].end(); ++mit)
      result[row].push_back(mit->second);
  }

  //however, result has dimensions [N][4]
  return result;

}



//convert a horz_matrix parse tree into an IMATRIX
//create an [N][4] matrix
template <typename TreeT>
IMATRIX horz_matrix(TreeT const& tree){
  
  IMATRIX result;
  std::map<char, std::vector<int> > result_map;
  std::map<char, std::vector<int> >::iterator mit;

  int matrix_length = tree.children[0].children.size() - 1;
  result.resize(matrix_length);

  for (size_t row = 0; row != tree.children.size(); ++row){
    TreeT const& rownode = tree.children[row];
    char nuc = *(rownode.children[0].value.begin());
    result_map[nuc] = std::vector<int>(matrix_length);
    
    for (int col = 0; col != matrix_length; ++col){
      std::string txt(rownode.children[col+1].value.begin(),
                      rownode.children[col+1].value.end());
      result_map[nuc][col] = atoi(txt.c_str());
    }
  }

  //this step iterates through the map in nucleotide key order (alphabetical)
  //thus ordering the rows of the results 
  
  for (int col = 0; col != matrix_length; ++col)
    for (mit = result_map.begin(); mit != result_map.end(); ++mit)
      result[col].push_back(mit->second[col]);

  return result;

}

template <typename IterT>
void assign_bool::operator()(IterT const& first,
                             IterT const& last) const {
  bool_val = std::string(first, last) == "true" ? true : false;
}


//assumes a flat tree structure containing name, data, and other things...
template <typename TreeT>
cis::pattern * MakePattern(TreeT const& tree){

  //iterate through the children of 'tree',
  //to get as many fields as possible
  std::string cluster_group;
  cis::dna_strand strand = cis::NON_STRANDED;
  int min_score = QUAL_SCORE;
  int total_hits = 0;
  std::string pattern_kind;
  IMATRIX matrix;
  map<std::string, cis::dnastring> consensus_pats;
  std::string name;
  std::map<string, cis::dnastring> consensi;
  typedef typename TreeT::const_tree_iterator tree_it;
    
  for (tree_it ch = tree.children.begin(); ch != tree.children.end(); ++ch){

    std::string txt(ch->value.begin(), ch->value.end());
    switch(ch->value.id().to_long()){
    case cluster_groupID: cluster_group = txt; break;
    case min_scoreID: min_score = atoi(txt.c_str()); break;
    case total_hitsID: total_hits = atoi(txt.c_str()); break;
    case strandID: strand = (cis::dna_strand)cis_enum::from_string(txt); break;
    case matrix_entryID: pattern_kind = "matrix"; name = txt; break;
    case vert_matrixID: matrix = vert_matrix(*ch); break;
    case horz_matrixID: matrix = horz_matrix(*ch); break;
    case collection_entryID: pattern_kind = "collection"; name = txt; break;
    case consensus_lineID: AddToConsensusMap(*ch, consensus_pats); break;
        
    default: {
      fprintf(stderr, "MakePattern: cannot understand this ID\n");
      exit(37);
    }; break;
      
    }
  }     
    
  hit_info hinfo(name, strand, 0);

  cis::pattern * pat = 
    new cis::pattern(NULL, NULL, hinfo, NULL, pattern_kind, 
                     cluster_group, 0, total_hits, min_score);
    
  if (pattern_kind == "matrix") {
    expandConsensusMatrix(matrix);

    //adjust the min-score so that it corresponds to the zero-adjustment
    pat->min_score -= zeroAdjustMatrix(matrix);
    pat->init_matrix(matrix);
  }

  else if (pattern_kind == "collection") pat->init_collection(consensus_pats);
  else {}

  return pat;
}



template<typename TreeT, typename IterT>
void make_pattern::operator()(TreeT const& tree,
                              IterT const& /*first*/,
                              IterT const& /*last*/) const {
  this->patterns.push_back(MakePattern<TreeT>(tree));
}




//parse the consensus data part of the map
template <typename TreeT>
void AddToConsensusMap(TreeT const& tree,
                       std::map<string, cis::dnastring> & consensus_map){

  typedef typename TreeT::const_tree_iterator tree_it;

  //tree_it name_node, cons_node;
  TreeT const& name_node = tree.children[0];
  TreeT const& cons_node = tree.children[1];
  consensus_map.insert
    (std::make_pair
     (std::string(name_node.value.begin(), name_node.value.end()),
      cis::ToDNAString(std::string(cons_node.value.begin(), 
                                   cons_node.value.end()))));
}



template <typename TreeT>
cis::cluster * MakeCluster(TreeT const& tree){
  
  string name, kind;
  int window = 0;
  int min_reqs = 0;
  diversity_cluster::TALLY reqs_list; //for diversity_cluster
  typedef typename TreeT::const_tree_iterator tree_it;

  for (tree_it ch = tree.children.begin(); 
       ch != tree.children.end(); ++ch){

    std::string txt(ch->value.begin(), ch->value.end());

    switch(ch->value.id().to_long()){
    case cluster_nameID: name = txt; break;
    case cluster_typeID: kind = txt; break;
    case cluster_windowID: window = atoi(txt.c_str()); break;
    case cluster_min_requirementsID: min_reqs = atoi(txt.c_str()); break;
    case cluster_requirementID: { 
        std::string clustname(ch->children[0].value.begin(),
                              ch->children[0].value.end());
        std::string sclustreq(ch->children[1].value.begin(),
                              ch->children[1].value.end());
        int clustreq = atoi(sclustreq.c_str());
       reqs_list[clustname] = clustreq;
    }; break;

//     case cluster_requirementsID: {
//       for (tree_it req = ch->children.begin();
//            req != ch->children.end(); ++req){
//       }
//     }; break; //handled by cluster_requirementID
      
    default: {
      fprintf(stderr, "MakeCluster: cannot understand this ID\n");
      //exit(37);
    }; break;
      
    }
  }     
  
  cis::cluster * ret = NULL;
  if (kind == "diversity"){
    ret = new diversity_cluster(name, window, reqs_list, min_reqs);
  } else if (kind == "specific") {
    fprintf(stderr, "Specific clustering not implemented yet...");
    exit(38);
  }

  return ret;
}


template <typename TreeT, typename IterT>
void make_cluster::operator()(TreeT const& tree,
                              IterT const& /*first*/, 
                              IterT const& /*last*/) const {
  clusters.push_back(MakeCluster<TreeT>(tree));
}

#endif // _PARSE_FUNCS_HPP
