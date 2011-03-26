#ifndef _MATCH_TRIE
#define _MATCH_TRIE

#include <map>
#include <set>
#include "dna.h"
#include "nuctrie.h"

using std::map;
using std::set;
using cis::NUC;


template <typename TAG>
class match_trie : public nuc_trie<TAG> {

	map<NUC, match_trie<TAG> *> children;
	TAG * _label;
	NUC _childbits;

 public:

	typedef bool (*MATCH_FCN)(NUC, NUC);
	typedef int (*LOSS_FCN)(NUC, NUC);

	static MATCH_FCN match;
	static LOSS_FCN loss;
	static int min_score;

	match_trie();
	~match_trie();

	TAG * label();
	NUC childbits();
	match_trie<TAG> * child(NUC);
	void insert(NUC const*const pat, TAG const & label);

	void match_all(NUC const*const pat, int, std::set<TAG *> & labels);
	void score_all(NUC const*const pat, int, std::set<TAG *> & labels, int cur_mismatch);
/* 	static void set_match_fcn(MATCH_FCN); */
/* 	static void set_loss_fcn(LOSS_FCN); */
/* 	static MATCH_FCN get_match_fcn(); */
/* 	static LOSS_FCN get_loss_fcn(); */
};


template <typename TAG>
bool (*match_trie<TAG>::match)(NUC query, NUC target);

template <typename TAG>
int (*match_trie<TAG>::loss)(NUC query, NUC target);

template <typename TAG>
int match_trie<TAG>::min_score;


//template <typename TAG> bool (match_trie<TAG>::* match)(NUC query, NUC target);

template<typename TAG>
match_trie<TAG>::match_trie(){
	label = NULL;
}


template<typename TAG>
match_trie<TAG>::~match_trie(){
	
	typename map<NUC, match_trie<TAG> *>::iterator cit;
	for (cit = children.begin(); cit != children.end(); cit++){
		delete (*cit).second;
	}
	children.clear();

	if (label != NULL) {
		delete label;
		label = NULL;
	}

}



template <typename TAG>
void match_trie<TAG>::match_all(NUC const*const target,
								int tsize,
								set<TAG *> & labels){

	typename map<NUC, match_trie<TAG> *>::iterator cit;
	if (label != NULL) labels.insert(label);

	if (tsize == 0) return;
	
	for (cit = children.begin(); cit != children.end(); cit++){
		if (match(cit->first, *target)){
			(*cit).second->match_all(target+1, tsize-1, labels);
		}
	}
}


//this assumes a 'score' field in the TAG, which is set temporarily
//this is in accordance with the 'early stopping criterion' whereby the
//scoring function is 
template <typename TAG>
void match_trie<TAG>::score_all(NUC const*const target, 
								int tsize,
								set<TAG *> & labels,
								int cur_score){

	typename map<NUC, match_trie<TAG> *>::iterator cit;
	if (label != NULL){
		label->score = cur_score; //this is only valid for one top-level complete call to score_all
		labels.insert(label);
	}

	if (tsize == 0) return;

	int next_score = 0;
	for (cit = children.begin(); cit != children.end(); cit++){
		next_score = cur_score + match_trie<TAG>::loss(cit->first, *target);
		if (next_score >= match_trie<TAG>::min_score)
			(*cit).second->score_all(target+1, tsize-1, labels, next_score);
	}
}


#endif //_MATCH_TRIE
