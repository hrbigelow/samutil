#ifndef _PWM_H
#define _PWM_H

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include "nuctrie.h"
#include "dna.h"

typedef nuc_trie<hit_info> hit_trie;

/* The idea is that each node will have a pwm_value associated with it,
   which is the matrix cell of the original PWM.
   Then, whenever it is traversed, its pwm_value is added to the 'score'
   variable.  Also, a minimum score is kept, and when the
 */

typedef std::vector<std::vector<int> > IMATRIX;

class pwm_trie : public hit_trie {

 public:
	typedef std::map<NUC, hit_trie *> MAP;
	typedef MAP::iterator IT;
	typedef MAP::value_type DATA;
	
	typedef pair<int, pwm_trie *> CELL; //score associated with the node...
	typedef std::vector<CELL> COL;
	typedef std::vector<COL> MAT;


	IMATRIX const* pwm_weights_;
	COL * _kids;

	hit_info const* label;

	int column;
	int ncolumns;
	int curscore;
	int minscore;
	int maxscore;

	pwm_trie(IMATRIX const* pwm_weights, COL * k, 
			 hit_info const* l, int c, int n, int s, int m) 
		: pwm_weights_(pwm_weights), _kids(k), label(l), 
		column(c), ncolumns(n), curscore(s), minscore(m) { }
	
	pwm_trie() 
		: pwm_weights_(NULL), _kids(NULL), label(NULL), 
		column(-1), ncolumns(0), curscore(0), minscore(0) { }
	
    ~pwm_trie();

	inline bool leaf() const { 
      assert((label != NULL) == (column == ncolumns - 1)); 
      return column == ncolumns - 1; 
    }

	int score(cis::dnastring const& hit_sequence) const;
	MAP children(NUC);
	void collect(std::set<hit_trie *> &, NUC const*const, int);
	inline void reset() { };
	inline int depth() const { return column; }
    void set_scores(int minscore, int maxscore);
    void set_curscore(int curscore);
	inline hit_info const* get_label() const { return label; }

 private:
    pwm_trie(pwm_trie const&);
	pwm_trie & operator=(pwm_trie const& f);

};


//problem:  no way to store the score, since the nodes are not one-to-one
//mapped with actual prefixes...note though that the scores *are* correct
//at the moment when they are stored...



//IMATRIX MakeMatrix(VNT const& dat, cis::dna_strand strand);
IMATRIX ReverseComplement(IMATRIX const& mat);

pwm_trie * MakePWM(IMATRIX const* pwm_weights, 
				   hit_info const* pat_label,
				   pwm_trie::MAT *);		



void PrintMatrix(IMATRIX const& matrix);
void PrintPWMMatrix(pwm_trie::MAT const& matrix);
void expandConsensusMatrix(IMATRIX & mat);
int zeroAdjustMatrix(IMATRIX & mat);

	



/* I thought I'd already solved this...
In any case, with a PWM, you don't actually insert any nodes at the beginning.  You merely
set a minimum score.  This will remain constant at the root node and will be passed on to
all created child nodes.

Children are created any time the children() function is called.  It will create on-demand all children
that score at or above the minimum score.  A second variable, positions_left, will be inherited by children,
but at one less than the parent.

The final part of the API will be that, when children at the zero level are created, their label field is initialized.



 */


#endif // _PWM_H
