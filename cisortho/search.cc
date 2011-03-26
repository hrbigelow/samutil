#include <string>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "search.h"
#include "nuctrie.h"
#include "pattern.h"
#include "index_trie_scan.h"
#include "hit.h"
#include "litestream.h"


//scans scantrie for pattern, returning the intervals found
LOCI Search(index_trie_scan<hit_trie> & scantrie,
            hit_trie * patnode,
            litestream & trie_stream,
            litestream & gpos_stream,
            int max_depth)
{

    LOCI ivals, ivals_filt;

    trie_stream.seekg(0, litestream::POS_BEGIN);
    gpos_stream.seekg(0, litestream::POS_BEGIN);

    //printf("flag state %i\n", trie_stream.rdstate());

    node_info disk_root_node = readNode(trie_stream);

    std::cout<<"Scanning...";
    scantrie.scan(trie_stream, disk_root_node, 0, patnode, 0, ivals);
    printf("found %i matching partial patterns.\n", (int)ivals.size());

    trie_stream.seekg(0, litestream::POS_BEGIN);
    node_info node = readNode(trie_stream);
    //ind_stream.close();

    printf("Filtering...");
    ivals_filt = scantrie.filter_by_suffix(ivals, gpos_stream, max_depth);
    printf("found %"PRId64" matching complete patterns.\n", ivals_filt.size());

    //   trie_stream.seekg(trie_position, litestream::POS_BEGIN);
    //   gpos_stream.seekg(gpos_position, litestream::POS_BEGIN);

    return ivals_filt;

}


int64_t TotalSize(index_trie_scan<hit_trie> & scantrie,
                  litestream & trie_stream )
{

    //std::ifstream ind_stream(scantrie.tree_file.c_str(), ios::in | ios::binary);

    node_info root = readNode(trie_stream);

    trie_stream.seekg(-root.size, litestream::POS_CURRENT);

    //ind_stream.close();

    return root.suffix_count;
}


int CountHits(LOCI const& loci)
{

    int count = 0;
    for (LOCI::const_iterator it = loci.begin(); it != loci.end(); it++)
    count += (*it).second.index_count;
    return count;
}


//returns the pair <score, num_extra_hits> given 
void FindThresholdScore(LOCI const& loci, int nhits, 
                        std::pair<int, int> *threshold)
{

    LOCI::const_iterator it;
    std::map<int, int> tally;
    std::map<int, int>::iterator mit;
    hit_trie * node;

    for (it = loci.begin(); it != loci.end(); ++it)
    {
        node = (*it).first;
        scan_interval const& sc = (*it).second;
        int nscore = - node->score(sc.sequence);
        tally[nscore]+= sc.index_count;
    }

    threshold->second = 0;
    //this is in numerically increasing order by key
    for (mit = tally.begin(); mit != tally.end(); ++mit)
    {
        if (threshold->second > nhits) break;
        threshold->first = - (*mit).first;
        threshold->second += (*mit).second;
    }
}



//makes new hits from ivals by scanning the position file,
//selecting only those hits residing on <dna>.
//max_depth tells us what the longest pattern found was...

//Should we assume that ivals only contains loci matching 'pat'?

//here we need the dna for each range of g coordinates...
//the most natural structure would be std::vector<pair<cis::dna_t const*, std::vector<int64_t> >

cis::REG_MAP makeHits(LOCI const& ivals,
                      cis::pattern *pat,
                      litestream & gpos_stream,
                      cis::dna_collection const& dnac,
                      cis::dna_collection::const_iterator start_dna,
                      cis::dna_collection::const_iterator end_dna,
                      int max_depth)
{
	
    cis::REG_MAP hits;

    //FILE * pos_stream = fopen(pos_file.c_str(), "r");
    //std::ifstream pos_stream(pos_file.c_str(), ios::in | ios::binary);

    hit_trie * patnode;

    cis::dna_t const* dna;

    for (LOCI::const_iterator iit = ivals.begin(); iit != ivals.end(); iit++)
    {

        patnode = (*iit).first;
        scan_interval const& scan = (*iit).second;

        std::vector<int64_t> gcoords =
        getGCoords(scan, gpos_stream, dnac, start_dna, end_dna);

        std::map<cis::dna_t const*, std::vector<int64_t> > loci = GetGroupedLoci(gcoords, dnac);

        //this is inconsistent with the idea of having only
        //terminal nodes at this point...
        hit_info const& terminal_label = *patnode->get_label();
        string hitname = 
        terminal_label.id.empty() ? 
        pat->label.id : 
        pat->label.id + '.' + terminal_label.id;

        cis::dna_strand hitstrand = 
        terminal_label.strand == cis::NON_STRANDED ? 
        pat->label.strand : 
        terminal_label.strand;
        
        int hitsize = terminal_label.size == 0 ? 
        pat->max_length : 
        terminal_label.size;

        for (std::map<cis::dna_t const*, std::vector<int64_t> >::iterator mit = loci.begin();
             mit != loci.end(); ++mit)
        {
            dna = mit->first;
            std::vector<int64_t>::iterator locus_iter;
            for (locus_iter = mit->second.begin(); locus_iter != mit->second.end();
                 ++locus_iter)
            {
                
                int64_t locus = (*locus_iter);
                cis::hit * h = new cis::hit(*dna, scan.sequence, locus, locus + hitsize,
                                            hitstrand, patnode->score(scan.sequence), -1, 
                                            hitname, &pat->cluster_group);
                hits[dna].insert(h);
            }
        }
    }
    
    //pos_stream.close();
    //fclose(pos_stream);

    return hits;
}
