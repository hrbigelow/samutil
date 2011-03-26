#include "index_trie_scan.h"
#include "nuctrie.h"

#include <cstdio>
#include <cassert>
#include <algorithm>


/*****************************************
index_trie_scan
******************************************/

std::vector<int64_t> getGCoords(scan_interval const& iv, 
                             litestream& pos,
                             cis::dna_collection const& dnac,
                             cis::dna_collection::const_iterator start_iterator,
                             cis::dna_collection::const_iterator end_iterator){
  
  int64_t gmin = dnac.ToIndex(start_iterator, 0);
  int64_t gmax = dnac.ToIndex(end_iterator, 0);

  std::vector<int64_t> gcoords;

  int64_t ibytes = iv.index_count * POS_BYTES;
  char * gpos = new char[ibytes];

  assert(pos.good());
  //fseeko(pos, 0, SEEK_END);
  //int64_t e = ftello(pos);
    
  //fseeko(pos, iv.index_start * POS_BYTES, SEEK_SET);
  //fread(gpos, 1, ibytes, pos);

  //pos.seekg();
  pos.seekg(iv.index_start * POS_BYTES, litestream::POS_BEGIN);
  pos.read(gpos, ibytes); //will this do the right thing?

  int64_t g;
  for (int64_t i=0; i < ibytes; i+=POS_BYTES){
    g = readInt(gpos+i, POS_BYTES);
    if (g >= gmin && g <= gmax) gcoords.push_back(g);
  }

  delete gpos;
  return gcoords;
}


std::vector<int64_t> getGCoords(scan_interval const& iv, litestream& pos){

  assert(pos.good());

  std::vector<int64_t> gcoords(iv.index_count);

  int64_t ibytes = iv.index_count * POS_BYTES;
  char * gpos = new char[ibytes];

  pos.seekg(iv.index_start * POS_BYTES, litestream::POS_BEGIN);
  pos.read(gpos, ibytes); //will this do the right thing?

  for (int64_t i=0; i < iv.index_count; i++) 
    gcoords[i] = readInt(gpos + (POS_BYTES*i), POS_BYTES);

  //sort them for optimized retrieval
  std::sort(gcoords.begin(), gcoords.end());

  delete gpos;
  return gcoords;
}


//gets the set of chromosomal start positions in a scan_interval
//that reside on 'dna'

//get the loci from a set of genomic coordinates
//the piece of dna is assumed known
std::map<cis::dna_t const*, std::vector<int64_t> > 
GetGroupedLoci(std::vector<int64_t> & gcoords, 
               cis::dna_collection const& dnac){
	
  std::map<cis::dna_t const*, std::vector<int64_t> > loci;

  std::pair<cis::dna_t const*, int64_t> locus;
  for (std::vector<int64_t>::const_iterator git = gcoords.begin(); 
       git != gcoords.end(); ++git){
    locus = dnac.FromIndex(*git);
    loci[locus.first].push_back(locus.second);
  }

  return loci;
}


// //gets loci spanning multiple pieces of DNA, presumed
// //to reside within the dna collection
std::vector<std::pair<cis::dna_t const*, int64_t> > 
GetIndividualLoci(std::vector<int64_t> const& gcoords,
                  cis::dna_collection const& dnac){
	 
  std::vector<std::pair<cis::dna_t const*, int64_t> > loci(gcoords.size());	
	 
  for (int i=0; i < (int)gcoords.size(); i++)
    loci[i] = dnac.FromIndex(gcoords[i]);
	 
  return loci;
}


