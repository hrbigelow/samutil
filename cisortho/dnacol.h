#ifndef _DNACOL_H
#define _DNACOL_H


#include <set>
#include <map>
#include <stdint.h>

#include "dna_types.h"

namespace std { template <typename T1, typename T2> class pair; }
namespace cis {
  
  class dna_t;
  struct less_dna_ptr;

  typedef std::set<cis::dna_t const*, less_dna_ptr> DNAS;

//class to provide a coordinate system and access to a set of pieces of dna, either chromosomes or contigs
 class dna_collection : public DNAS {

	 int64_t total_length_;
	 std::map<cis::dna_t const*, int64_t> offsets;
	 std::map<int64_t, cis::dna_t const*> offsets_inv;
	 typedef std::map<cis::dna_t const*, int64_t>::const_iterator OIT;
	 typedef std::map<int64_t, cis::dna_t const*>::const_iterator OIIT;

 public:
	 int64_t ToIndex(cis::dna_t const*, int64_t) const; //returns a unique index for the position and piece of DNA
	 int64_t ToIndex(DNAS::const_iterator, int64_t) const; //returns a unique index for the position and piece of DNA
	 std::pair<cis::dna_t const*, int64_t> FromIndex(int64_t) const; //returns the position and piece of DNA given the unique index
	 void add(std::string const species, 
              std::string const& fasta_dir,
              std::string const& fasta_file);
	 int64_t num_bases() const;
     void open_dnas();
	 void calc_offsets();

     //initialize and make the dna_collection ready for retrieving and scanning
     void make_ready(char const* dna_index_file, std::string dna_directory);

	 dna_collection(DNAS const& d = DNAS());


	 //calculate the total length of dnas in the range [start_dna, end_dna)
	 int64_t TotalLength(DNAS::const_iterator dna_start, 
					  DNAS::const_iterator dna_end) const;

 };

 cis::dna_t const* GetDNAByName(DNAS const& dnas, std::string const& species, std::string const& sdna);


}


#endif // _DNACOL_H
