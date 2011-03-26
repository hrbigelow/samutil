#include <stdint.h>

#include "dnastats.h"
#include "dna_scanning.h"
#include "dna_types.h"

using namespace std;

namespace cis { class dna_t; }
//is called at every window position during the initial indexing scan for IterateDNAWindows
//if the prefix of pat matches prefix, it is inserted.  otherwise it's not
class collectSyms : public scanStore {

public:
	//dnastring prefix;
	bool found[16];
	
	void operator()(void *root, cis::NUC const*const pat, int patsize, int64_t /*gpos*/, cis::dna_t const* /*dna*/){
		set<cis::NUC> & nucs = *((set<cis::NUC> *)root);
		for (int i=0; i < 16; i++) found[i] = false;
		for (int i=0; i < patsize; i++){ found[(int)pat[i]] = true; }
		for (int i=0; i < 16; i++) if (found[i]){ nucs.insert(cis::NUC(i)); }
	}

	collectSyms() { for (unsigned int i=0; i < 16; i++) found[i] = false; }
};



class collectNmers : public scanStore {

public:
	
	void operator()(void *root, cis::NUC const*const pat, int patsize, int64_t /*gpos*/, cis::dna_t const* /*dna*/){
		set<cis::dnastring> & nmers = *((set<cis::dnastring> *)root);
		nmers.insert(cis::dnastring(pat, patsize));
	}

	collectNmers() { }
};



set<cis::NUC> dnastats::Alphabet(cis::dna_t const& d, int chunk_size){
	set<cis::NUC> nucs;
	collectSyms scan;
	IterateDNAWindows(&scan, &nucs, chunk_size, chunk_size, chunk_size, &d);
	return nucs;
}



set<cis::dnastring> dnastats::Spectrum(cis::dna_t const& d, int nmer_size, int chunk_size){
	set<cis::dnastring> nmers;
	collectNmers scan;
	IterateDNAWindows(&scan, &nmers, nmer_size, nmer_size, chunk_size, &d);
	return nmers;
}
