#ifndef ICISJ_GRAMMAR_HPP_
#define ICISJ_GRAMMAR_HPP_

//#define BOOST_SPIRIT_DEBUG
//#define BOOST_SPIRIT_DEBUG_FLAGS BOOST_SPIRIT_DEBUG_FLAGS_TREES
#define BOOST_SPIRIT_DEBUG_PRINT_SOME 60

#include <climits>
#include <map>

//#include <boost/spirit.hpp>

#include <boost/spirit/core.hpp>
#include <boost/spirit/tree/parse_tree.hpp>
#include <boost/spirit/tree/ast.hpp>
#include <boost/spirit/utility/loops.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/utility/chset.hpp>
#include <boost/spirit/utility/lists.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/actor/assign_actor.hpp>

using namespace boost::spirit;

typedef char const*         iterator_t;
typedef tree_match<iterator_t> parse_tree_match_t;
typedef parse_tree_match_t::tree_iterator iter_t;


enum nodeIDS {

  inputID, cluster_groupID, total_hitsID, min_scoreID,
  strandID, pattern_kindID, pattern_nameID, vert_matrixID,
  matrix_entryID, collection_entryID, horz_matrixID, pattern_propID, collectionID,
  clusterID, cluster_nameID, cluster_typeID, cluster_windowID,
  consensus_lineID, cluster_min_requirementsID, cluster_requirementsID,
  cluster_requirementID

};

  



template <typename OutputT>
struct cis_grammar : public grammar<cis_grammar<OutputT> > {

  cis_grammar(OutputT & output_) : output(output_) {}
 

  template <typename ScannerT>
  struct definition {


    definition(cis_grammar const& g) {

      colon          =  discard_node_d [ ch_p(':') ];
      alnum_word     =  token_node_d[lexeme_d[+alnum_p]];
      graph_word     =  token_node_d[lexeme_d[+graph_p]];
      file_path      =  token_node_d[lexeme_d[+(alnum_p | chset_p("./-_"))]];
      

      //Global settings ********************************************
      global_settings = 
        data_directory | merge_overlap | species
        | position_file | tree_file | dna_index_file;
        
        
      species = str_p("species") >> colon >>
        token_node_d[lexeme_d[+graph_p]][push_back_a(g.output.species)];

//       species_list = str_p("species") >> colon >> 
//         list_p(token_node_d[lexeme_d[+(graph_p - ',')]]
//                [push_back_a(g.output.species)], ch_p(','));
      
      data_directory = str_p("data_directory") >> colon >>
        graph_word[assign_a(g.output.data_directory)];

      merge_overlap = str_p("merge_overlapping_clusters") >> colon >>
        as_lower_d [ str_p("true") | str_p("false") ]
        [g.output.set_merge_overlap];
      
      position_file = str_p("position_file") >> colon >>
        graph_word[assign_a(g.output.posfile)];

      tree_file = str_p("tree_file") >> colon >>
        graph_word[assign_a(g.output.treefile)];

      dna_index_file = str_p("dna_index_file") >> colon >>
        graph_word[assign_a(g.output.dnaindexfile)];



      // global limits on input space **********************************
      global_limits = max_total_hits | max_collection_size 
        | max_pwm_length | min_pwm_length | min_collection_score | max_patterns;

        max_total_hits = str_p("max_total_hits") >> colon >> int_p
        [assign_a(g.output.max_total_hits)];

      max_collection_size = str_p("max_collection_size") >> colon >> int_p
        [assign_a(g.output.max_collection_size)];

      max_pwm_length = str_p("max_pwm_length") >> colon >> int_p
        [assign_a(g.output.max_pwm_length)];

      min_pwm_length = str_p("min_pwm_length") >> colon >> int_p
        [assign_a(g.output.min_pwm_length)];

      min_collection_score = str_p("min_collection_score") >> colon >> int_p
        [assign_a(g.output.min_collection_score)];

      max_patterns = str_p("max_patterns") >> colon >> int_p
        [assign_a(g.output.max_patterns)];

      
    
    
      //Pattern for finding hits ****************************************
      pattern = 
        (matrix_entry >> vert_matrix >> +pattern_prop)
        | (matrix_entry >> horz_matrix >> +pattern_prop)
        | (collection_entry >> +consensus_line >> +pattern_prop)
        //| discard_node_d[alignment_tag >> colon] >> pattern_name >> alignment
        ;
 
      matrix_entry = discard_node_d[str_p("matrix") >> colon] >> 
        token_node_d[lexeme_d[+graph_p]];
      
      
     //Vertical Matrix
      vert_matrix_header = repeat_p(4)[nucleotide_sym];
      vert_matrix_line   = repeat_p(4)[int_p];
      vert_matrix_lines  = +vert_matrix_line;
      vert_matrix        = vert_matrix_header >> vert_matrix_lines;
    
      //Horizontal Matrix
      horz_matrix_line   = nucleotide_sym >> discard_node_d[colon] >> +int_p;
      horz_matrix        = repeat_p(4)[horz_matrix_line];

      pattern_name       = graph_word;
    
      consensus_line   = graph_word >> consensus_seq;

      collection_entry = discard_node_d[str_p("collection") >> colon] >>
        token_node_d[lexeme_d[+graph_p]];

      nucleotide_sym =  chset_p("ACGT");
      consensus_sym  =  chset_p("ACMGRSVTWYHKDBN");

      consensus_seq  =  token_node_d[lexeme_d[+chset_p("ACMGRSVTWYHKDBN")]];
      
      //       exact_nmer     =  token_node_d[ lexeme_d[ +nucleotide_sym ] ];
      //       alignment        =  +(exact_nmer);
    
      //Pattern properties *****************************************
      pattern_prop  = cluster_group | total_hits | min_score | strand;
      cluster_group = discard_node_d [ str_p("cluster_group") >> colon ] >> 
        leaf_node_d[lexeme_d[+graph_p]];

      //pattern_name;

      total_hits    = discard_node_d [ str_p("total_hits") >> colon ] >> int_p;
      min_score     = discard_node_d [ str_p("min_score") >> colon ] >> int_p;
      strand        = discard_node_d [ str_p("strand") >> colon ] >> 
        (str_p("POS") | str_p("NEG"));
    
    
      //Clustering criteria **********************************************
      cluster = +(cluster_name | cluster_type | cluster_window |
                  cluster_min_requirements | cluster_requirement);

      cluster_name  = discard_node_d [ str_p("cluster") >> colon] >> 
        leaf_node_d[lexeme_d[+graph_p]];

      cluster_type  = discard_node_d [str_p("cluster_type") >> colon] >> 
        (str_p("diversity") | str_p("specific"));

      cluster_window = discard_node_d[str_p("cluster_window") >> colon] >> int_p;
      cluster_min_requirements = 
        discard_node_d [str_p("cluster_min_requirements") >> colon] >> int_p;

//       cluster_requirements =
//         discard_node_d[str_p("cluster_requirements_list") >> colon] >> 
//         //gen_pt_node_d[list_p(cluster_requirement, ch_p(','))];
//       list_p(cluster_requirement, discard_node_d[ch_p(',')]);

      cluster_requirement = 
        discard_node_d[str_p("cluster_requirement") >> colon] >>
        leaf_node_d[lexeme_d[+graph_p]] >> leaf_node_d[int_p];


    
    
      //The Main Input rule--just a large collection of top-level rules here...
      input =
        +(global_limits 
          | access_node_d[pattern][g.output.pattern_maker]
          | global_settings
          | access_node_d[cluster][g.output.cluster_maker]
          );


//       BOOST_SPIRIT_DEBUG_NODE(vert_matrix_header);
//       BOOST_SPIRIT_DEBUG_NODE(vert_matrix_line);
//       BOOST_SPIRIT_DEBUG_NODE(vert_matrix_lines);
//       BOOST_SPIRIT_DEBUG_NODE(vert_matrix);
//       BOOST_SPIRIT_DEBUG_NODE(horz_matrix_line);
//       BOOST_SPIRIT_DEBUG_NODE(horz_matrix);
      BOOST_SPIRIT_DEBUG_TRACE_NODE_NAME(pattern, "pattern", BOOST_SPIRIT_DEBUG_TRACENODE);
//       BOOST_SPIRIT_DEBUG_NODE(pattern_name);
//       BOOST_SPIRIT_DEBUG_NODE(pattern_prop);
//       BOOST_SPIRIT_DEBUG_NODE(cluster_group);
//       BOOST_SPIRIT_DEBUG_NODE(total_hits);
//       BOOST_SPIRIT_DEBUG_NODE(min_score);
//       BOOST_SPIRIT_DEBUG_NODE(strand);
//       BOOST_SPIRIT_DEBUG_NODE(matrix_entry);
      
//       BOOST_SPIRIT_DEBUG_NODE(alnum_word);
//       BOOST_SPIRIT_DEBUG_NODE(graph_word);
//       BOOST_SPIRIT_DEBUG_NODE(file_path);
//       BOOST_SPIRIT_DEBUG_NODE(colon);
//       BOOST_SPIRIT_DEBUG_NODE(species_list);
//       BOOST_SPIRIT_DEBUG_NODE(collection);
//       BOOST_SPIRIT_DEBUG_NODE(collection_entry);
//       BOOST_SPIRIT_DEBUG_NODE(consensus_line);
//       BOOST_SPIRIT_DEBUG_NODE(consensus_seq);

      BOOST_SPIRIT_DEBUG_NODE(cluster);
//       BOOST_SPIRIT_DEBUG_NODE(cluster_window);
//       BOOST_SPIRIT_DEBUG_NODE(cluster_min_requirements);
//       BOOST_SPIRIT_DEBUG_NODE(cluster_requirements);
//       BOOST_SPIRIT_DEBUG_NODE(cluster_requirement);

      
      
    };
  
    rule<ScannerT> const& start() const { return input; }

    rule<ScannerT> input;

    rule<ScannerT> colon, graph_word, file_path, global_settings, 
      merge_overlap, data_directory,
      global_limits, max_total_hits, max_collection_size, max_pwm_length, 
      min_pwm_length, min_collection_score, species, max_patterns, pattern, 
      vert_matrix_header, vert_matrix_line, vert_matrix_lines, horz_matrix_line,
      position_file, tree_file, dna_index_file;

    rule<ScannerT> alnum_word;


    rule<ScannerT, parser_context<>, parser_tag<matrix_entryID> > matrix_entry;
    rule<ScannerT, parser_context<>, parser_tag<collection_entryID> > collection_entry;
    
    rule<ScannerT, parser_context<>, parser_tag<vert_matrixID> > vert_matrix;
    rule<ScannerT, parser_context<>, parser_tag<horz_matrixID> > horz_matrix;
    rule<ScannerT, parser_context<>, parser_tag<collectionID> > collection;

    rule<ScannerT, parser_context<>, parser_tag<pattern_nameID> > pattern_name;

    rule<ScannerT, parser_context<>, parser_tag<consensus_lineID> > consensus_line;
  
    rule<ScannerT> pattern_prop;
    rule<ScannerT, parser_context<>, parser_tag<cluster_groupID> > cluster_group;
    rule<ScannerT, parser_context<>, parser_tag<total_hitsID> > total_hits;
    rule<ScannerT, parser_context<>, parser_tag<min_scoreID> > min_score;
    rule<ScannerT, parser_context<>, parser_tag<strandID> > strand;

    rule<ScannerT, parser_context<>, parser_tag<pattern_kindID> > 
    collection_tag, alignment_tag;

    //    rule<ScannerT, parser_context<>, parser_tag<alignmentID> > alignment;
  
    rule<ScannerT> nucleotide_sym, consensus_sym, consensus_seq;
  

    //Cluster rules
    rule<ScannerT, parser_context<>, parser_tag<clusterID> > cluster;
    rule<ScannerT, parser_context<>, parser_tag<cluster_nameID> > cluster_name;
    rule<ScannerT, parser_context<>, parser_tag<cluster_typeID> > cluster_type;
    rule<ScannerT, parser_context<>, parser_tag<cluster_windowID> > cluster_window;
    rule<ScannerT, parser_context<>, parser_tag<cluster_min_requirementsID> > cluster_min_requirements;
    rule<ScannerT, parser_context<>, parser_tag<cluster_requirementsID> > cluster_requirements;
    rule<ScannerT, parser_context<>, parser_tag<cluster_requirementID> > cluster_requirement;

  };

  OutputT & output;
  
};

#endif

