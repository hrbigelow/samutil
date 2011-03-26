#include <cassert>
#include <vector>
#include <algorithm>

#include "pwm.h"

pwm_trie::~pwm_trie(){
  for (size_t i = 0; i < (*_kids).size(); ++i){
    if (i) _kids->clear();
    else delete (*_kids)[i].second;
  }
}


//return the children such that the resulting
//prefix still passes the score threshold.
//The policy will be to return an empty set if we're at the end

/* Return all nodes (usually 1 or zero) that match 'n' and
   also pass the cumulative score threshold.
   The _children variable caches all qualifying children that
   are seen up to this point.  They are produced on-demand as
   the pattern is seen.

   The returned value is a subset of these which match the nucleotide...
 
*/
int pwm_trie::score(cis::dnastring const& hit_sequence) const {
  int score = 0;

  assert(hit_sequence.size() >= (*pwm_weights_).size());

  for (int i=0; i < (int)pwm_weights_->size(); ++i)
    score += (*pwm_weights_)[i][hit_sequence[i]];
  return score;
}


//For whatever child that is returned, the curscore
//of that child is correctly set, assuming the curscore
//of the calling node is set correctly...
pwm_trie::MAP pwm_trie::children(NUC nuc){

  if (leaf()) return pwm_trie::MAP();
  pwm_trie::MAP ret;
  int nextscore = curscore + (*_kids)[static_cast<int>(nuc)].first;
  if (nextscore >= minscore){
    pwm_trie * child = (*_kids)[static_cast<int>(nuc)].second;
    child->curscore = nextscore;
    if (child->leaf() && curscore > maxscore){ }
    else { ret[nuc] = child; }
  }
  return ret;
}



//we need some form of collect that collects non-terminal nodes...
//something that's re-entrant;  i.e. something that you can
void pwm_trie::collect(set<nuc_trie<hit_info> *> & pats, 
					   NUC const*const target, int tsize){

  if (column == ncolumns-1) pats.insert(this);
  else {
    MAP kids = children(*target);
    for (IT ch = kids.begin(); ch != kids.end(); ++ch)
      (*ch).second->collect(pats, target+1, 0);
  }	
}


//recursively set the minimum score of all nodes
void pwm_trie::set_scores(int minscore, int maxscore){
  this->minscore = minscore;
  this->maxscore = maxscore;
  for (size_t i = 0; i < (*_kids).size(); ++i){
    if (i) {
      (*_kids)[i].second->minscore = minscore;
      (*_kids)[i].second->maxscore = maxscore;
    } else {
      (*_kids)[i].second->set_scores(minscore, maxscore);
    }
  }
}


void pwm_trie::set_curscore(int curscore){
  curscore = curscore;
}


typedef std::vector<int> ICOLUMN;
typedef std::vector<ICOLUMN> IMATRIX;


//TODO: add label initialization routine here.
pwm_trie * MakePWM(IMATRIX const* pwm_weights, 
				   hit_info const* pat_label, 
				   pwm_trie::MAT *const grid){

  int ncols = pwm_weights->size();
  pwm_trie::MAT & gr = *grid;

  gr.resize(ncols+1); //allocates empty vectors
  gr[ncols] = pwm_trie::COL(); //empty last column...

  hit_info const* current_label;
	
  for (int c = ncols - 1; c >= 0; c--){
    gr[c].resize(16);
    current_label = (c == ncols - 1) ? pat_label : NULL;
    assert((*pwm_weights)[c].size() == 16);
    for (int n = 0; n < 16; n++)
      gr[c][n] = 
        pwm_trie::CELL((*pwm_weights)[c][n], 
                       new pwm_trie(pwm_weights, & gr[c+1],
                                    current_label, c, ncols, 0, 0));
  }

  return new pwm_trie(pwm_weights, & gr[0], NULL, -1, ncols, 0, 0);
}

//maps the four (A,C,G,T) scores to the 16 consensus
//nucleotides.  Takes the minimum (worst) score of each
//For example, if the column were 5,2,1,4, then [AC] would be 2, N woudl be 1
//This ensures that patterns such as NNNNNN don't score high, but that
//patterns containing only one or 2 Ns are still permitted


void expandConsensusMatrix(IMATRIX & mat){

  int v[4];
  for (size_t i=0; i < mat.size(); ++i){
    std::vector<int> & r = mat[i];
    std::vector<int> r2(16, MIN_SCORE);
    for (int j=0; j < 16; ++j){
      v[0] = j & cis::A ? r[0] : MAX_SCORE;
      v[1] = j & cis::C ? r[1] : MAX_SCORE;
      v[2] = j & cis::G ? r[2] : MAX_SCORE;
      v[3] = j & cis::T ? r[3] : MAX_SCORE;
      r2[j] = std::min(v[0], std::min(v[1], std::min(v[2], v[3])));
      if (r2[j] == MAX_SCORE) r2[j] = MIN_SCORE;
    }
    mat[i] = r2;
  }
}


//add a constant to all matrix columns such that their max is zero
int zeroAdjustMatrix(IMATRIX & mat){

  //int maxv = MIN_SCORE;
  int total_decrement = 0;
  for (unsigned int i=0; i < mat.size(); i++){
    std::vector<int> & r = mat[i];
    int maxv = *std::max_element(r.begin(), r.end());
    for (unsigned int j=0; j < r.size(); j++) r[j] -= maxv;
    total_decrement += maxv;
  }
  return total_decrement;

}


IMATRIX ReverseComplement(IMATRIX const& mat){

  IMATRIX rc(mat);
  reverse(rc.begin(), rc.end());
  IMATRIX::iterator mit;
  for (mit = rc.begin(); mit != rc.end(); mit++){
    ICOLUMN const& col = (*mit);
    ICOLUMN rccol(col.size());
    for (unsigned int i=0; i<col.size(); i++)
      rccol[i] = col[(int)cis::CompNUC(i)];
    (*mit) = rccol;
  }
  return rc;
}



// IMATRIX MakeMatrix(VNT const& dat, cis::dna_strand strand){

//   IMATRIX mat = jf::ParseMatrix(jf::getv(dat,"data"), &jf::exti);
//   expandConsensusMatrix(mat);
//   zeroAdjustMatrix(mat);
//   if (strand == cis::NEG) mat = ReverseComplement(mat);
	
//   return mat;
	
// }



void PrintMatrix(IMATRIX const& matrix){
  for (size_t c = 0; c != matrix.size(); ++c){
    for (size_t r = 0; r != matrix[c].size(); ++r){
      if (r == 0) printf("%i", matrix[c][r]);
      else printf(", %i", matrix[c][r]);
    }
    printf("\n");
  }
}


void PrintPWMMatrix(pwm_trie::MAT const& matrix){
  for (size_t c = 0; c != matrix.size(); ++c){
    for (size_t r = 0; r != matrix[c].size(); ++r){
      if (r == 0) printf("%i", matrix[c][r].first);
      else printf(", %i", matrix[c][r].first);
    }
    printf("\n");
  }
}
