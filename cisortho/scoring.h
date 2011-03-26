#ifndef _SCORING_H
#define _SCORING_H

#include "enum.h"
#include "dna.h"
#include "match_trie.h"
#include "boost/bind.hpp"
#include "boost/function.hpp"
#include <utility>
#include <cmath>
#include <vector>
#include <set>
#include <string>
#include <map>

using std::pair;
using std::vector;
using std::set;
using std::string;
using std::map;

namespace cis 
{
//typedef boost::function<pair<int, int> (NUC const*)> SCORING;
	typedef vector<vector<int> > PWM;
	typedef vector<vector<float> > PFM;
	
	
	bool iupac2nuc_match(NUC q, NUC t);
	int iupac2nuc_loss(NUC q, NUC t);




class scoring {
	static unsigned int ctr;

 public:
	const unsigned int id;
	const string name;
	unsigned int size;
	enum scoretype { EXACT, TRIE_MATCH, TRIE_SCORE, MATRIX, SET, CONSENSUS } stype;
	virtual pair<int, int> scores(NUC const*) = 0;

	bool exact;
	//pair<int, int> scores(const char * seq) { return func(seq); }
	scoring(string const& n, const unsigned int s, scoretype st) 
		: id(ctr++), name(n), size(s), stype(st) {}
	scoring() : id((unsigned int)HUGE_VAL), name(""), size(0) {}
	//scoring operator=(scoring const& sc) { assert(0); }
	virtual ~scoring() {}

};


class consensus_score : public scoring {

	dnastring cons;
	dnastring cons_rc;
	unsigned int max_mismatch;

 public:
	consensus_score(string const&, char const*, unsigned int = 0);
	pair<int, int> scores(NUC const* n);
	
};


class consensus_set_score : public scoring {

	set<dnastring> cons;

 public:
	consensus_set_score(string const&, set<dnastring>);
	pair<int, int> scores(NUC const* n);
	
};


class pwm_score : public scoring {

	int max_possible_score;
	PWM pwm;
	PWM rev;

 public:
	
	pwm_score(string const&, PWM const&);
	pair<int, int> scores(NUC const* n);
};


class pfm_score : public scoring {

	PFM pfm;
	PFM rev;

 public:
	
	pfm_score(string const&, PFM const&);
	pair<int, int> scores(NUC const* n);
};


class set_score : public scoring {

	set<dnastring> forward;
	set<dnastring> revcomp;

 public:

	set_score(string const&, set<dnastring>);
	pair<int, int> scores(NUC const* n);
};


 //typedef pair<string, dna_strand> SENSE_ID;
typedef map<hit_info, dnastring> SENSE_PATS;
 typedef match_trie<hit_info> TRIE_T;

 class trie_match : public scoring {

	 TRIE_T * pattern_trie;
	 TRIE_T::MATCH_FCN match_fcn;

 public:
	 
	 set<hit_info *> ids_found;
	 
	 trie_match(string const& name, SENSE_PATS pats, TRIE_T::MATCH_FCN);
	 ~trie_match();
	 pair<int, int> scores(NUC const*);
	 
 };
 
 
 class trie_score : public scoring {

	 TRIE_T * pattern_trie;

	 const TRIE_T::LOSS_FCN loss_fcn;
	 const int min_score;
	 const int initial_score;
	 
 public:
	 
	 set<hit_info *> ids_found;
	 trie_score(string const& name, SENSE_PATS pats, TRIE_T::LOSS_FCN, int, int);
	 ~trie_score();

	 pair<int, int> scores(NUC const*);

 };



template <typename V> vector<vector<V> > 
RevComp(vector<vector<V> > & mat){
	vector<vector<V> > rev(mat.size()); //recursive copy?
	int pN = mat.size() - 1;
	int pC;
	for (unsigned int c = 0; c < mat.size(); c++){
		vector<int> &pos = mat[c];
		pC = pN - c;
		vector<int> posrev(pos.size());
		posrev[0] = pos[3];
		posrev[1] = pos[2];
		posrev[2] = pos[1];
		posrev[3] = pos[0];
		posrev[4] = pos[4];
		rev[pC] = posrev;
	}
	return rev;
}


dnastring RevComp(dnastring const&);
string RevComp(string const& s);
void TranslateSeq(NUC *const iseq, char const*const seq, unsigned int const size);
void Complement(NUC *const icomp, NUC *const iseq, unsigned int const size);
string FromDNAString(dnastring const& d);
dnastring ToDNAString(string const& s);
int descending(const void *,const void *);

scoring * LoadScoringFcn(const char *fpwm, char const*, const char *indir);
PWM LoadGenericPWM(bfs::path const&) throw (string&);
PWM LoadHMMERPWM(bfs::path const&) throw (string&);


} //namespace cis


#endif // _SCORING_H
