#ifndef _DNA_SCANNING_H
#define _DNA_SCANNING_H

#include <stdint.h>

#include "dna_types.h"

namespace cis { class dna_t; };

class scanStore {
public:
	virtual void operator()(void *, cis::NUC const*const, 
                            int, int64_t, cis::dna_t const*) = 0;
	virtual ~scanStore(){}
};


void IterateDNAWindows(scanStore * fcn, void * store, int minwinsize, 
					   int maxwinsize, int chunk_size, cis::dna_t const* dna);

#endif // _DNA_SCANNING_H
