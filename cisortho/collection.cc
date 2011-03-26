#include "collection.h"

using cis::dnastring;

typedef collection_trie CT;

CT::~collection_trie(){

  MAP::iterator it;
  for (it = _children.begin(); it != _children.end(); ++it) delete (*it).second;
  if (label != NULL) {
    delete label;
    label = NULL;
  }
}

//score policy for this will be the following:
//the score of the current node is assumed correct.
//the scores of the child nodes are updated to be in the
//correct relation to the current node, or are untouched
//if they don't qualify.
//therefore, you will never see scores that are below the
//minimum score...


//create and return the subset of children whose mismatches
//fall below the number of maximum mismatches
//we want the base of this collection trie to be a superset
//of the nucleotide
CT::MAP CT::children(NUC nuc){
  MAP sel;
  CT * ch;
  for (IT chit = _children.begin(); chit != _children.end(); ++chit){
    ch = (CT *)(*chit).second;
    ch->mismatch = (nuc & ch->_base) && !(nuc & ~ch->_base) ? 
      mismatch : mismatch+1;
    if (ch->mismatch <= max_mismatch) sel.insert(*chit);
  }
  return sel;
}



void CT::insert(NUC const*const pat, int psize, hit_info const& lid){

  if (psize == 0) { 
    set_label(new hit_info(lid));
    return;
  }
  IT it = _children.find(*pat);
  CT * child;
  if (it == _children.end()){
    child = (CT *) (_children[*pat] = 
                    new collection_trie(*pat, 0, max_mismatch, _depth+1, NULL));
  }
  else child = (CT *)it->second;
	
  return child->insert(pat+1, psize-1, lid);
}


void CT::insert(dnastring const& d, hit_info const & label){
  return insert(d.c_str(), d.size(), label);
}



//stores all patterns in subtree of <this> that match target with
//up to 'max_mismatch' mismatches.
void CT::collect(set<nuc_trie<hit_info> *> & pats,
				 NUC const*const target,
				 int tsize){
	
  if (mismatch > max_mismatch) return;
	
  if (label != NULL) pats.insert(this);
	
  if (tsize == 0) return;

  CT * ch;
  for (IT it = _children.begin(); it != _children.end(); it++){
    ch = (CT *)(*it).second;
    ch->mismatch = mismatch;
    if (! (*target & (*it).first)) ch->mismatch++;
    ch->collect(pats, target+1, tsize-1);
  }


}

	
int CT::max_depth(int cd) const {

  CT * ch;
  cd = std::max(depth(), cd);
  for (MAP::const_iterator it = _children.begin(); it != _children.end(); it++){
    ch = (CT *)(*it).second;
    cd = std::max(ch->max_depth(cd), cd);
  }
  return cd;
	
}	




collection_trie * MakeColTree(map<string, dnastring> const& pats, int max_mismatch){
  collection_trie * t = new collection_trie(cis::X, max_mismatch, 0, 0, NULL);
  map<hit_info, dnastring> sites;
  map<hit_info, dnastring>::iterator it;
  pair<string, dnastring> sd;
  dnastring dpat, dpat_rc;
  map<string, dnastring>::const_iterator pats_iter;
  for (pats_iter = pats.begin(); pats_iter != pats.end(); ++pats_iter){
    dpat = (*pats_iter).second;
    dpat_rc = cis::RevCompDNA(dpat);
    if (dpat == dpat_rc){
      sites[hit_info((*pats_iter).first, cis::POS_NEG, dpat.size())] = dpat;
    } else {
      sites[hit_info((*pats_iter).first, cis::POS, dpat.size())] = dpat;
      sites[hit_info((*pats_iter).first, cis::NEG, dpat.size())] = dpat_rc;
    }				
  }			
  for (it = sites.begin(); it != sites.end(); it++)
    t->insert((*it).second, (*it).first);
  return t;
}
