#include <set>

#include "dna_scanning.h"
#include "dna_types.h"

namespace cis { 
  class dna_t; 
}

namespace dnastats {
  std::set<cis::NUC> Alphabet(cis::dna_t const& d, int chunk_size);
  std::set<cis::dnastring> Spectrum(cis::dna_t const& d, int nmer_size, int chunk_size);
}
