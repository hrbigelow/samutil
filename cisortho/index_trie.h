#ifndef _INDEX_TRIE_H
#define _INDEX_TRIE_H

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <stdint.h>

#include "dna_types.h"

namespace cis { class dna_collection; }

class litestream;

const int POS_BYTES = 5;

struct less_nuc {
	bool operator()(cis::NUC const& n1, cis::NUC const& n2) const {
		return n1 < n2;
	}
};


template <typename T>
bool strnsame(T const*s1, T const* s2, int n)
{
  do { if (*s1++ != *s2++) return false;
  } while (--n != 0);
  return true;
}


class index_node_label;

class index_trie {

 public:

	enum node_type { leaf = 0, uleaf, branch, stub };

	int64_t nbytes;
	int64_t nsuffixes;
	bool is_merge_node;

	typedef std::map<cis::NUC, index_trie *, less_nuc> CMAP;
	typedef CMAP::value_type CHILD;
	typedef CMAP::iterator ITER;

	static unsigned int max_suffix_load;
	static void set_max_suffix_load(unsigned int i) { max_suffix_load = i; }

	CMAP children;
	index_node_label *label;
	index_trie();
	~index_trie();
	void build_subtree(cis::dnastring, cis::dna_collection const&, int, int, int);
	void build_toptree(std::string const& treedir, std::string const& posdir);
	void insert(int64_t position, cis::NUC const*const pat, short size, short left);
	void insert(int64_t, cis::dnastring, short);

	void insertfile(cis::NUC const*const, int,
                    std::string const& treefile, int64_t,
                    std::string const& posfile, int64_t);
	void count();
	void push();
	char header(cis::NUC);
	void write(cis::NUC, std::ostream &, std::ostream &);
	node_type type();
	bool empty();
			
};

const unsigned short int SUFFIX_BYTES = 63;
const unsigned short int MAX_NODE_BYTES = 64;

union NUCpair {
	struct {
		cis::NUC N1 : 4;
		cis::NUC N2 : 4;
	};
	char NN;
};


union Header {

	struct {
		unsigned int S : 3;
		unsigned int T : 3;
		index_trie::node_type N: 2;
	};
	char dat;
};


class index_node_label {
 public:
	//typedef std::set<pair<int, cis::dnastring> >::iterator SIT;
	typedef std::pair<int64_t, cis::dnastring> SPOS;
	typedef std::vector<SPOS> SFX;
	typedef std::vector<SPOS>::iterator SIT;
	std::string treefile;
	int64_t treesize;
	std::string posfile;
	int64_t poscount; //the number of loci
	SFX suffixes;
	index_node_label() : treefile(""), treesize(0), posfile(""), poscount(0) {}
};



struct node_info {

	int size;
	index_trie::node_type type;
	cis::NUC nuc;
	int64_t suffix_count;
	int64_t subtree_size;
	NUCpair suffix[SUFFIX_BYTES];
	int suffix_bytes;
	int child_count;

	node_info(Header, char const*, unsigned int);

	friend std::ostream& operator<<(std::ostream&, node_info const&);

};


node_info readNode(litestream &);
int64_t readInt(char const*const ptr, int nbytes);


std::ostream& printBytes(char const*, std::ostream&, litestream&, int);


#endif // _INDEX_TRIE_H
