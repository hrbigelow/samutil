/************************
Routines for producing patterns from parsed data.
*************************/
#include <vector>
#include <numeric>

#include "pattern.h"
#include "collection.h"
#include "pwm.h"
#include "dna.h"
#include "enum.h"

using cis::dnastring;

namespace {

  hit_info MakeHitInfo(string name, string pat){

	dnastring dpat = cis::ToDNAString(pat);
	dnastring dpat_rc = cis::RevCompDNA(dpat);
	cis::dna_strand ds = dpat == dpat_rc ? cis::POS_NEG : cis::POS;
	return hit_info(name, ds, pat.size());

  }

//   dnastring MakeDNAString(VNT v){
//     return cis::ToDNAString(jf::exts(v));
//   }

} // namespace


namespace cis
{
	
  pattern::pattern(pattern const& p) :
    ptr(NULL),
    ptr2(NULL),
    label(p.label),
    t(NULL),
    kind(p.kind),
    cluster_group(p.cluster_group),
    max_length(p.max_length),
    total_hits(p.total_hits),
    min_score(p.min_score)
  { }
	
  
  pattern::pattern(void * ptr,
                   void * ptr2,
                   hit_info label,
                   nuc_trie<hit_info> * ntrie,
                   std::string const& kind,
                   std::string const& cg,
                   int max_length,
                   int total_hits,
                   int min_score) :
    ptr(ptr), ptr2(ptr2), label(label), t(ntrie),
    kind(kind), cluster_group(cg), max_length(max_length),
    total_hits(total_hits), min_score(min_score){}
                  


  //make a copy of 'matrix', which is assumed to be
  //a [N][16] PWM matrix over the 16 consensus symbols.
  //it is assumed to be zero-adjusted (zero is the best score achievable)
  void pattern::init_matrix(IMATRIX const& matrix){

    pwm_trie::MAT * grid = new pwm_trie::MAT();
    IMATRIX * pwm_weights = new IMATRIX(matrix);

    ptr = grid;
    ptr2 = pwm_weights;
    label.size = pwm_weights->size();
    pwm_trie *pt = MakePWM(pwm_weights, & label, grid);
    pt->set_scores(min_score, 0);
    max_length = pwm_weights->size();
    t = pt;
  }
  
  
  void pattern::init_collection(std::map<std::string, cis::dnastring> const& 
                                consensus_pats){
    
    collection_trie *ct = MakeColTree(consensus_pats, -min_score);
    t = ct;
    max_length = ct->max_depth(0);
  }
  
  
  pattern::~pattern(){
    if (t != NULL){ delete t; t = NULL; }
    if (ptr != NULL) { 
      if (kind == "matrix"){ 
        delete static_cast<IMATRIX *>(ptr); 
        delete static_cast<IMATRIX *>(ptr2);
        ptr = NULL; 
        ptr2 = NULL;
      }
      else {
        cerr<<"Don't know what ptr points to!"<<endl;
        exit(53);
      }
    }
  }


  //make a new reverse-complemented pattern from the ready-to-use, provided
  //pattern 'pat'
  pattern * MakeRCPattern(pattern const* pat){
    if (pat->kind != "matrix"){
      fprintf(stderr, "Only matrix patterns can be reverse complemented");
      exit(36);
    }
    pattern *pat_rc  = new pattern(*pat); 

    //preserves 'NON_STRANDED' and 'POSNEG' settings...
    if (pat->label.strand == cis::POS) pat_rc->label.strand = cis::NEG;
    else if (pat->label.strand == cis::NEG) pat_rc->label.strand = cis::POS;
    else pat_rc->label.strand = pat->label.strand;

    IMATRIX pwm_weights = 
      ReverseComplement(* static_cast<pwm_trie *>(pat->t)->pwm_weights_);

    //this will redundantly call expandConsensusMatrix and 
    //zeroAdjustMatrix, but will do no harm...(yes it will...)
    pat_rc->init_matrix(pwm_weights);
    return pat_rc;
  }


} // namespace cis
