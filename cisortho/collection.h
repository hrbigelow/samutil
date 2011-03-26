#include <map>
//#include "dna_types.h"
#include "nuctrie.h"

class collection_trie : public nuc_trie<hit_info> {

	typedef std::map<cis::NUC, nuc_trie<hit_info> *> MAP;
	typedef MAP::iterator IT;
	typedef MAP::value_type DATA;

	cis::NUC _base;
	MAP _children;
	int max_mismatch;
	int mismatch;
	int _depth;
	hit_info * label;

 public:

	collection_trie(cis::NUC n, int mm, int m, unsigned int d, hit_info * l) : 
		_base(n), max_mismatch(mm), mismatch(m), _depth(d), label(l) { }

	~collection_trie();
	inline bool leaf() const { return _children.empty(); }
	inline int score(cis::dnastring const& /*hit_sequence*/) 
		const { return max_mismatch - mismatch; }

	MAP children(cis::NUC n);
	void insert(cis::NUC const*const pat, int size, hit_info const & label);
	void insert(cis::dnastring const& d, hit_info const & label);
	void collect(set<nuc_trie<hit_info> *> &, cis::NUC const*const, int);
	inline void reset() { mismatch = 0; }
	inline int depth() const { return _depth; }
	inline void set_label(hit_info * l){ label = l; }
	inline hit_info * get_label() const { return label; }

	int max_depth(int) const;
	
};


collection_trie * MakeColTree(std::map<std::string, cis::dnastring> const&, int);
