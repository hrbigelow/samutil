#ifndef _INDEX_TRIE_SCAN_H
#define _INDEX_TRIE_SCAN_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <cstring>
#include <sys/stat.h>

#include "dna.h"
#include "dnacol.h"
#include "index_trie.h"
#include "litestream.h"



//describes the set of loci in the trie having 'sequence' and residing at
//various positions.
//this object will be paired with a 'pat node', which itself has an intrinsic
//position relative to the patterns here. 
struct scan_interval {
  int64_t index_start; //what is this?  this is the cumulative number of suffixes passed up so far
  int64_t index_count; //the number of suffixes in the current subtrie
  cis::dnastring sequence;
  scan_interval(int64_t s, int64_t c, cis::NUC * seqarray, int seqsize) : 
    index_start(s), index_count(c), sequence(seqarray, seqsize) {}
};


template <typename TPAT>
class index_trie_scan {

 public:
  cis::dna_collection dnac;
  std::string tree_file;
  std::string pos_file;

  int max_depth;

  cis::NUC *prefix;

  index_trie_scan(cis::dna_collection const& d,
                  std::string const& t,
                  std::string const& p,
                  int m) : 
    dnac(d), tree_file(t), pos_file(p), max_depth(m) {
    prefix = new cis::NUC[max_depth+1];
    struct stat64 stat_buf;
    if (stat64(tree_file.c_str(), &stat_buf) != 0){
      cerr<<"index trie file "<<tree_file<<" doesn't exist"<<endl;
      exit(50);
    }
    if (stat64(pos_file.c_str(), &stat_buf) != 0){
      cerr<<"locus file "<<pos_file<<" doesn't exist"<<endl;
      exit(50);
    }

  }

  index_trie_scan() : 
    tree_file(std::string()), 
    pos_file(std::string()), max_depth(0), prefix(NULL) {
  }
	
  ~index_trie_scan(){ 
    if (prefix != NULL) { delete prefix; prefix = NULL; }
  }
	
  void scan(litestream & index, 
            node_info const&,
            int, //suffix index
            TPAT *, 
            int, //prefix size
            std::vector<std::pair<TPAT *, scan_interval> >&);

  std::vector<std::pair<TPAT *, scan_interval> >
    filter_by_suffix(std::vector<std::pair<TPAT *, scan_interval> >&, 
                     litestream & gpos_stream, int);
	
};


std::vector<int64_t> getGCoords(scan_interval const& iv, 
                                  //FILE * pos,
                                  litestream& pos, 
                                  cis::dna_collection const&,
                                  cis::dna_collection::const_iterator start_iterator,
                                  cis::dna_collection::const_iterator end_iterator);

std::vector<int64_t> getGCoords(scan_interval const& iv, litestream& pos);
std::map<cis::dna_t const*, std::vector<int64_t> > GetGroupedLoci(std::vector<int64_t> & gcoords, 
											  cis::dna_collection const&);
std::vector<std::pair<cis::dna_t const*, int64_t> > GetIndividualLoci(std::vector<int64_t> const& gcoords,
												  cis::dna_collection const& dnac);


//stores in ivals the means for generating the complete sequence(s) that matches TPAT.
//the prefix is necessary because patnode doesn't store it, furthermore it may be degenerate.
//the suffix can either be retrieved directly from the disknode, or scanned from the DNA
//if there are multiple suffixes or if the single suffix is too short...

//should we encapsulate the rest of the search in 'store'?
//that way.  in a sense, 

template <typename TPAT>
void store(std::vector<std::pair<TPAT *, scan_interval> >& ivals,
		   TPAT * patnode,
		   node_info const& dnode,
		   int64_t index_start,
		   cis::NUC const*const pat,
		   int patsize){

  //make the pat and the patnode compressed suffix into one pattern
  //and store it as a dnastring
  const int seqsize = patsize + dnode.suffix_bytes * 2;
  cis::NUC * sequence = new cis::NUC[seqsize];
  memcpy(sequence, pat, patsize);
  cis::NUC * sequence_extension = sequence + patsize;
  for (int i=0; i < dnode.suffix_bytes; i++){
    sequence_extension[2*i] = dnode.suffix[i].N1;
    sequence_extension[2*i+1] = dnode.suffix[i].N2;
  }

  ivals.push_back(std::make_pair(patnode, 
                            scan_interval(index_start, 
                                          dnode.suffix_count,
                                          sequence, seqsize)));
  delete sequence;
}



//find patterns represented by patnodes in the index tree
//the <patnodes> can represent
//this must represent 
template <typename TPAT>
void index_trie_scan<TPAT>::scan(litestream & ind_stream,
                                 node_info const& disk_parent_node,
                                 int current_prefix_num,
                                 TPAT * search_parent_node,
                                 int prefix_sequence_index,
                                 std::vector<std::pair<TPAT *, scan_interval> >& ivals){
	
  //node_info disk_child_node; //disk child
  TPAT *search_child_node; //pat child

  std::map<cis::NUC, TPAT *> search_child_nodes;
  typename std::map<cis::NUC, TPAT *>::iterator nit;
  int child_count = disk_parent_node.child_count;
	
  while (child_count--){
		
    node_info disk_child_node = readNode(ind_stream); //get the disk child

    if (disk_child_node.type != index_trie::stub){
      search_child_nodes = search_parent_node->children(disk_child_node.nuc);
      if (search_child_nodes.empty()){
        ind_stream.seekg(disk_child_node.subtree_size - disk_child_node.size, 
                         litestream::POS_CURRENT);
        current_prefix_num += disk_child_node.suffix_count;
        continue;
      }
    }
    prefix[prefix_sequence_index] = disk_child_node.nuc;

    switch(disk_child_node.type){
    case index_trie::leaf: 
    case index_trie::uleaf:
      for (nit = search_child_nodes.begin(); 
           nit != search_child_nodes.end(); nit++)
        store(ivals, (*nit).second, disk_child_node, current_prefix_num,
              prefix, prefix_sequence_index+1);
      current_prefix_num += disk_child_node.suffix_count;
      break;

    case index_trie::stub:
      if (search_parent_node->leaf())
        store(ivals, search_parent_node, disk_parent_node,
              current_prefix_num, prefix, prefix_sequence_index);
      current_prefix_num += disk_child_node.suffix_count;
      break;

    case index_trie::branch: {
      int64_t pos = ind_stream.tellg(); //just after reading the disk node header
      for (nit = search_child_nodes.begin(); 
           nit != search_child_nodes.end(); 
           ++nit){
        ind_stream.seekg(pos, litestream::POS_BEGIN);
        search_child_node = (*nit).second;
        if (search_child_node->leaf()) {
          cis::dnastring dprefix(prefix, prefix_sequence_index+1);
/*           printf("%i\t%s\n", (int)disk_child_node.suffix_count, */
/*                  cis::FromDNAString(dprefix).c_str()); */
          store(ivals, search_child_node, disk_child_node,
                current_prefix_num, prefix, prefix_sequence_index+1);

        } else {
          scan(ind_stream, disk_child_node, current_prefix_num, 
               search_child_node, prefix_sequence_index+1, ivals);
        }
      }
      //done with this entire subtree.
      ind_stream.seekg(pos+disk_child_node.subtree_size-disk_child_node.size,
                       litestream::POS_BEGIN);
      current_prefix_num += disk_child_node.suffix_count;
      break;
    }
    }
  }
}

/* Scan the DNA for non-terminal pattern nodes, obtaining a filtered
   interval list of complete intervals.  The final output list can be
   used to count the total number of nodes, and also to retrieve the
   loci information from the position file, and make hits.  For the
   non-terminal pattern nodes, the loci must be queried here too, and
   then the DNA file itself.
*/

template <typename TPAT>
std::vector<std::pair<TPAT *, scan_interval> > index_trie_scan<TPAT>::filter_by_suffix(
    std::vector<std::pair<TPAT *, scan_interval> >& ivals,
    litestream & gpos_stream,
    int max_depth){
	
  typename std::set<TPAT *>::iterator ist;

  std::vector<std::pair<TPAT *, scan_interval> > ret;

  TPAT * pat;
  cis::dna_t const* dna;
  int64_t chr_pos;
  char * seq_extension = new char[max_depth];
  cis::NUC * dseq = new cis::NUC[max_depth];

  int64_t gindex_start;

  typename std::vector<std::pair<TPAT *, scan_interval> >::iterator iit;

  for (iit = ivals.begin(); iit != ivals.end(); ++iit){

    pat = (*iit).first;
    scan_interval & ival = (*iit).second;
    int ssize = ival.sequence.size();

/*     std::vector<int64_t> gcoords = getGCoords(ival, gpos_stream); */
/*     std::vector<std::pair<cis::dna_t const*, int64_t> > loci = */
/*       GetIndividualLoci(gcoords, dnac); */

		
    //*****************************************
    //*****************************************

    if (pat->leaf()) {
      assert(pat->get_label() != NULL);
      ret.push_back(*iit);
    } else { //pat is not a terminal, search dna file and evaluate all subpatterns
			
      std::vector<int64_t> gcoords = getGCoords(ival, gpos_stream);
      std::vector<std::pair<cis::dna_t const*, int64_t> > loci =
        GetIndividualLoci(gcoords, dnac);

      gindex_start = ival.index_start;

      int i = 0;
      int extension_bytes;
      int direct_read_bytes;
      for (std::vector<std::pair<cis::dna_t const*, int64_t> >::iterator it = loci.begin(); 
           it != loci.end(); ++it, ++i){

        dna = (*it).first;
        chr_pos = (*it).second;
        if (chr_pos + max_depth > dna->length()) continue;

        //std::cout<<dna->name<<":"<<chr_pos<<endl;
        extension_bytes = std::max(max_depth - ssize, 0);
        direct_read_bytes = std::min(ssize, max_depth);
				
        memcpy(dseq, ival.sequence.c_str(), direct_read_bytes);

        if (extension_bytes){
          dna->source().seekg(dna->seek_start_pos+chr_pos+ssize, litestream::POS_BEGIN);
          dna->source().read(seq_extension, extension_bytes);
        }

        cis::TranslateSeq(dseq + ssize, seq_extension, extension_bytes);

        std::set<TPAT *> terminals;

        //a sufficient piece of information would be simply the complete DNA sequence,
        //plus an offset
        pat->collect(terminals, dseq + pat->depth(), max_depth - pat->depth());

        pat->reset();

        assert(terminals.size() <= 1);

        for (ist = terminals.begin(); ist != terminals.end(); ++ist){
          scan_interval
            terminal_interval(gindex_start + i, 1, dseq, ssize+extension_bytes);
          ret.push_back(std::make_pair(*ist, terminal_interval));
        }
				
      }
    }
  }

  delete seq_extension;
  delete dseq;
  //pos_stream.close();

  return ret;
}





#endif // _INDEX_TRIE_SCAN_H
