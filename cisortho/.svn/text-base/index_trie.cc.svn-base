#include <cmath>
#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstring>

#include "index_trie.h"
#include "dna.h"
#include "dna_scanning.h"
#include "dnacol.h"

#include "litestream.h"


//return the minimum number of bytes
//necessary to store the int
const int64_t bits255 = 255;

int64_t byteMask[] = { 
    bits255, 
    bits255<<8, 
    bits255<<16, 
    bits255<<24, 
    bits255<<32,
    bits255<<40,
    bits255<<48,
    bits255<<56
};

int intBytes(unsigned int v)
{
    int b;
    for (b = 3; b >= 0; b--) if (v & byteMask[b]) break;
    return b+1;
	
}

//returns the minimum number of bytes to store this long integer
int64_t longBytes(int64_t v)
{
    int64_t b;
    for (b = 7; b >= 0; b--) if (v & byteMask[b]) break;
    return b+1;
}

//this will be affected by little and big endian notation...
void writeInt(ostream & f, int64_t val, int sz)
{
    if (sz) f.write((char *)&val, sz);
}


//read and interpret an integer value using only nbytes of data
int64_t readInt(char const*const ptr, int nbytes)
{
    int64_t n = 0;
    memcpy((char *)&n, ptr, nbytes);
    return n;
}
	


index_trie::index_trie() : is_merge_node(false), label(new index_node_label)
{
}

index_trie::~index_trie()
{
    for (ITER cit = children.begin(); cit != children.end(); cit++)
    delete (*cit).second;
    delete label;
}


//cannot be called with n = 0

//is called at every window position during the initial indexing scan for IterateDNAWindows
//if the prefix of pat matches prefix, it is inserted.  otherwise it's not
class insertMatching : public scanStore 
{

public:

    cis::dna_collection dnac;
    int mindepth;
    int maxdepth; //only used to check
    cis::NUC *pfx;
    int plen;

    void operator()(void *root, cis::NUC const*const pat, int patsize, int64_t gpos, cis::dna_t const* dna)
    {
        if (strnsame(pfx, pat, plen))
        {
            ((index_trie *)root)->insert
            (dnac.ToIndex(dna, gpos), pat+plen, patsize-plen, mindepth);
        }
    }

    insertMatching(cis::dnastring p, cis::dna_collection dc, int md, int xd) : 
        dnac(dc), mindepth(md), maxdepth(xd)
    {
        pfx = new cis::NUC[p.size()];
        plen = p.size();

        memcpy(pfx, p.c_str(), plen);
        if (mindepth < plen)
        {
            cerr<<"Minimum tree depth ("<<mindepth<<") must be >= Prefix Length ("<<plen<<")"<<endl;
            exit(53);
        }
        if (maxdepth < mindepth)
        {
            cerr<<"Maximum tree depth must be >= minimum tree depth of "<<mindepth
                <<" given."<<endl;
            exit(53);
        }
    }

    ~insertMatching()
    {
        if (pfx != NULL)
        { 
            delete pfx; pfx = NULL; 
        }
    }
};


/*
//is called at every window position during the initial indexing scan for IterateDNAWindows
//if the prefix of pat matches prefix, it is inserted.  otherwise it's not
class insertMatching : public scanStore {

public:

dna_collection dnac;
int mindepth;
int maxdepth; //only used to check
cis::NUC *pfx;
int plen;
index_trie * root;
cis::NUC * pat;
int patsize;
int64_t dna_offset;
	

void operator()(int winpos)
{
if (strnsame(pfx, pat + gpos, plen) && 
(char *)memchr(pat, cis::X, mindepth) == NULL)
((index_trie *)root)->insert
(dna_offset + winpos, pat+plen, patsize-plen, mindepth);
}

insertMatching(dnastring p, dna_collection dc, int md, int xd) : 
dnac(dc), mindepth(md), maxdepth(xd)
{
pfx = new cis::NUC[p.size()];
plen = p.size();

strncpy(pfx, p.c_str(), plen);
if (mindepth < plen)
{
cerr<<"Minimum tree depth ("<<mindepth<<") must be >= Prefix Length ("<<plen<<")"<<endl;
exit(53);
}
if (maxdepth < mindepth)
{
cerr<<"Maximum tree depth must be >= minimum tree depth of "<<mindepth
<<" given."<<endl;
exit(53);
}
}

~insertMatching()
{
if (pfx != NULL)
{ delete pfx; pfx = NULL; }
}
};
*/


 //builds an index tree of all <N>mer patterns beginning with <prefix> from <dnac>
 //reading <chunk_size> pieces of dna at a time for scanning
void index_trie::build_subtree(cis::dnastring prefix, cis::dna_collection const& dnac,
							   int mindepth, int maxdepth, int chunk_size)
{
	
    //if ((char *)memchr(prefix.c_str(), cis::X, prefix.size()) != NULL) return;
    insertMatching scan_fcn(prefix, dnac, mindepth, maxdepth); //this will alert any problem
    cis::dna_t const* dna;
    for (cis::dna_collection::iterator dnac_iter = dnac.begin(); 
         dnac_iter != dnac.end(); ++dnac_iter)
    {
        dna = (*dnac_iter);
        IterateDNAWindows(&scan_fcn, (void *)this, mindepth, maxdepth, 
                          chunk_size, dna);
    }
}

bool validHeader(Header h)
{

    return true;
    // 	return
    // 		h.N == index_trie::leaf && ! h.T && h.S || //leaf must have zero subtree size but positive suffix size
    // 		h.N == index_trie::branch && h.T && h.S || //branch must have positive subtree and suffix size
    // 		h.N == index_trie::single_suffix && h.L || //single suffix has non-zero size
    // 		h.N == index_trie::multiple_suffix && ! h.T && h.S; //multiple leaf must have zero subtree size and positive suffix size
}


/*

  node information:

  leaf          00000SSS N SSSSSSS
  uleaf         01000000 N X 1234567... 0 to 255

  branch        10TTTSSS N C TTTTTTT SSSSSSS

  stub          11000SSS SSSSSSS

  A leaf may have N = cis::nonNUC.  This signals that the leaf contains
  the suffixes belonging to its parent node.

  The transformation from the 


  branch        11TTTSSS N C TTTTTTT SSSSSSS X 1234567... 0 to 255  This is unnecessary, probably?

  branch        11TTTSSS N C TTTTTTT SSSSSSS XXXX   (XXXXXXXX is the number of truncated suffixes.

  The suffixes that reside on a branch will be all zero, and they won't be leaves because they won't
  be long enough.

  There only need be 4 bytes worth, amply sufficient (4,000,000,000) to store the truncated loci
  count for a given suffix...

  If this is the case, we should push patterns down?

  But, then, we need to implement push...

  //how to signal that a branch contains one or more suffixes?

  because there is no way to predict the full size of the header from the


  merge         <the on-disk subtrie is written as raw bytes here>

  The leaf, branch and root nodes represent actual tree nodes
  The multi-suffix and single-suffix represent the suffix vector of the
  previous node having one or many suffixes.  They always appear immediately after
  the node header of the node containing them.

*/

//it would be good to have the node_info constructor not have to copy anything.
//That is, it should initialize its members during construction using scanning...

//construct node information from a header and the rest of the node

node_info::node_info(Header h, char const* n, unsigned int rest) : 
    size(rest+1), type(h.N), suffix_bytes(0), child_count(0)
{

    char const* np = n;
	
    switch (type)
    {
    case index_trie::leaf: 
        nuc = *np++;
        suffix_count = readInt(np, h.S); np  +=  h.S;
        subtree_size = size;
        break;

    case index_trie::uleaf:
        nuc = *np++;
        suffix_count = 1;
        subtree_size = size;
        suffix_bytes = readInt(np++, 1);
        memcpy((char *)suffix, np, suffix_bytes);
        break;

    case index_trie::branch:
        nuc = *np++;
        child_count = readInt(np, 1); np++;
        subtree_size = readInt(np, h.T); np += h.T;
        suffix_count = readInt(np, h.S); np += h.S;
        break;

    case index_trie::stub:
        nuc = cis::z;
        suffix_count = readInt(np, h.S); np += h.S;
        subtree_size = size;
        break;
    }
	
}

	
//calculates the size of the chunk of info, then reads that and converts it into the node info...
//reads and returns the node information
node_info readNode(litestream & is)
{

    Header H;
    char node[MAX_NODE_BYTES];

    is.get(H.dat);

    //first index is type, second is whether it has a suffix
    //order is leaf, uleaf, branch, stub
    //inner order is 'no suffix', 'suffix'
    static unsigned int node_sizes[] = { 1, 2, 2, 0 };
		
    unsigned int rest = node_sizes[H.N] + H.S + H.T;

    is.read(node, rest);

    if (H.N == index_trie::uleaf)
    {
        int x = readInt(node+rest-1, 1);
        is.read(node+rest, x);
        rest += x;
    }

    return node_info(H, node, rest);

}



ostream& operator<<(ostream &o, node_info const& n)
{

    char const* nt = NULL;
    switch(n.type)
    {
    case index_trie::leaf: nt = "Leaf"; break;
    case index_trie::uleaf: nt = "ULef"; break;
    case index_trie::branch: nt = "Brnc"; break;
    case index_trie::stub: nt = "Stub"; break;
    }
	
    o<<nt<<' '<<n.child_count;
    // 	for (int i=0; i < 16; i++)
    //{
        // 		char c = (ch & m<<i) ? nucs[i] : '-';
        // 		o<<c;
        // 	}

        // 	o<<'\t'
        // 	 <<n.subtree_size<<'\t'
        // 	 <<n.size<<'\t'
        // 	 <<n.suffix_count;

        // 	for (int i=0; i < n.suffix_bytes; i++)
    //{
            // 		if (i ==0) o<<'\t';
            // 		o<<nucs[(int)n.suffix[i].N1]
            // 		 <<nucs[(int)n.suffix[i].N2];
            // 	}

    return o;
}


ostream& printBytes(char const* zero_one, ostream &o, litestream & is, int sz)
{
    
    char * node = new char[sz];
    is.read(node, sz);
    for (int i=0; i < sz; i++)
    {
        if (i > 0) o<<' ';
        for (int b=7; b >= 0; b--)
        {
            char j = zero_one[node[i]>>b & 1];
            o<<j;
        }
    }
    delete node;
    return o;
}


unsigned int index_trie::max_suffix_load;

/*
  we want to initially insert a pattern until it runs out, so the leaf
  will have an empty suffix, or until it becomes a single leaf somewhere
  is expressed as a unique prefix plus a suffix stored in the leaf
  the 'left' variable, when non-zero, overrides the behavior and ensures that
  the pattern will continue on down the trie, creating nodes as needed, until
  at the prefix is at least 'left' long.
  Since there is no test for the number of suffixes at a leaf, there is here
  no way to ensure that multiple suffixes get extended...
  In fact, we don't necessarily want multiple leaves to be eliminated.  Just
  that more than one suffix sharing a prefix

*/

//add a pattern defined by pat, size to tree starting at root
//will only add it to depth necessary to produce a unique leaf,
//or if the pattern runs out, a leaf with multiple identical patterns

//inserts the pattern by laying down as 

//but how do we insure that there will only be multiple leaves at the
//level of the max?

/* insert:  inserts a pattern defined by <pat>,<size> into the index_trie.
   Guarantees:

   1.  the inserted pattern will be stored in the trie as a prefix +
   suffix.  The prefix will be defined by the path of nodes in the
   trie.  The suffix will be stored as a dnastring in the 'suffixes'
   field of the last node.

   2.  A newly inserted pattern will be inserted to a depth such that
   the suffix itself does not share any prefix with other suffixes at
   that node.

   2a. The pattern will be inserted at least to a depth of max(size, left)

   3.  The pattern is guaranteed to be inserted as far into the
   existing trie as possible without creating any new nodes.  That is,
   the suffix could not be pushed down any further without creating a
   new node.

   Gotchas: 4.  The insertion of a new pattern may invalidate the
   depth maximality of any previous suffixes.  It does so because,
   upon encountering a leaf node with one or more suffixes, the new
   suffix may share a prefix with one or more of them.  In that case,
   Guarantee 2 mandates that it avoid sharing that prefix, and be
   pushed one level down.  This will create a new node, and thus
   invalidate Guarantee 3 for the clashing pattern.

   5.  Any inserted pattern will only invalidate suffixes below the
   point of pattern insertion.  Therefore, once all patterns are
   inserted, their suffixes can be updated by re-inserting each suffix
   in top-down order.  This is because, when a given suffix is
   encountered, if it is valid, it will remain valid (because there
   will be no more insertion events above it).  If it is invalid, it
   will be re-inserted at that point, and cannot invalidate any
   suffixes at that level or higher.
*/


/*
  In general, what should be the criterion for when to keep pushing vs when to allow
  the accumulation of suffixes at a given node?
  Suppose, if there are more than two suffixes, push them, otherwise let it slide!


*/
void index_trie::insert(int64_t pos, cis::NUC const*const pattern, 
                        short pattern_size, short min_pattern_size)
{
	
    if (pattern_size == 0)
    label->suffixes.push_back(make_pair(pos, cis::dnastring()));

    else {
        cis::NUC nuc(*pattern);
        ITER it = children.find(nuc);
        if (it != children.end())
        {
            (*it).second->insert(pos, pattern+1, pattern_size-1, min_pattern_size-1);
        } else if (min_pattern_size > 0)
        {
            index_trie *leaf = children[nuc] = new index_trie();
            leaf->insert(pos, pattern+1, pattern_size-1, min_pattern_size-1);
        } else {
            if (label->suffixes.size() > index_trie::max_suffix_load)
            {
                index_trie *leaf = children[nuc] = new index_trie();
                leaf->insert(pos, pattern+1, pattern_size-1, min_pattern_size-1);
            } else {
                label->suffixes.push_back(make_pair(pos, 
                                                    cis::dnastring(pattern, pattern_size)));
            }
        }
    }
}


void index_trie::insert(int64_t pos, cis::dnastring pattern, short min_pattern_size)
{
    return insert(pos, pattern.c_str(), pattern.size(), min_pattern_size);
}


//similar to inserting a nucleotide pattern, this inserts the pattern
//and creates a MERGE node at the end, to be used for merging individual files into one.
void index_trie::insertfile(cis::NUC const*const pat, int size, 
                            string const& tf, int64_t tsz,
                            string const& pf, int64_t pct)
{

    if (size > 0)
    {
        cis::NUC nuc(*pat);
        ITER it = children.find(nuc);
        if (it == children.end())
        {
            index_trie *leaf = new index_trie();
            children[nuc] = leaf;
        }
        children[nuc]->insertfile(pat+1, size-1, tf, tsz, pf, pct);
    }
    else {
        label->treefile = tf;
        label->treesize = tsz;
        label->posfile = pf;
        label->poscount = pct;
        is_merge_node = true;
    }

}

//calculate how many bytes it will take to store both the number of bytes
//in question plus a header that stores that number of bytes


/*
  if nv == v then nv == I(x + nv) == I(x + I(x + nv)) = I(x + I(x + I(x + nv)))
*/

//calculates the quantity
//x + intBytes(x + intBytes(x + intBytes(x + ...)))
//If x is the fixed number of bytes in the subtree plus the node
int64_t calcBytes(int64_t x)
{

    int64_t v = longBytes(x);
    int64_t nv = longBytes(x + v);
    while (nv != v)
    { 
        v = nv; nv = longBytes(x + v);
    }
    return x + nv;
}
	

//this tallies the counts of branches, suffixes, single leaves, and multiple leaves
//in the tree, by returning the child and collecting the stats from it

//this counts both the actual number of suffixes in the nodes
//and their children cumulatively, but also
//predicts the number of bytes that each node (including it's
//contained subtree) will take up.
//this second part is tricky because the number of bytes actually
//depends on the sizes of the sub-numbers themselves, thus
//introducing a somewhat circular dependence...
void index_trie::count()
{

    if (is_merge_node)
    {
        nbytes = label->treesize;
        nsuffixes = label->poscount;
        return;
    }

    int64_t s = (int64_t)label->suffixes.size();

    switch(type())
    {

    case leaf: nsuffixes = s; nbytes = 2 + longBytes(nsuffixes); break;
    case uleaf: 
        nsuffixes = 1; 
        nbytes = 3 + label->suffixes[0].second.size() / 2;
        break;

    case branch: {

        nsuffixes = s;

        int64_t sub_bytes = 0;
        for (ITER ch = children.begin(); ch != children.end(); ch++)
        {
            (*ch).second->count();
            sub_bytes += (*ch).second->nbytes;
            nsuffixes += (*ch).second->nsuffixes;
        }

        //add the size of what will be the stub node
        if(s) sub_bytes += 1 + longBytes(s);

        //the tricky part
        int64_t fixed_bytes = 3 + longBytes(nsuffixes) + sub_bytes;
        nbytes = calcBytes(fixed_bytes); //this step adds the T field
        break;
    }

    case stub:
        cerr<<"Shouldn't have a stub node here."<<endl;
        break;
    }
}


/*
//recursively go down the tree, and push any nonempty suffixes
//that are on branches down until they are on leaves
void index_trie::push()
{

static pair<int64_t, cis::dnastring> del(-1, cis::dnastring());

//control iterator invalidation...
label->suffixes.reserve(label->suffixes.size() * 2);

//find the current end, after the call to reserve
index_node_label::SIT sit, end = label->suffixes.end();
for (sit = label->suffixes.begin(); sit != end; sit++)
{

if ((*sit).second.empty()) continue; //nothing to do...
pair<int64_t, cis::dnastring> ps = *sit; //copy
		
(*sit).first = -1;
(*sit).second.clear();
		
insert(ps.first, ps.second.c_str(), ps.second.size(), 0);
		
}

sit = remove(label->suffixes.begin(), label->suffixes.end(), del);
label->suffixes.erase(sit, label->suffixes.end());
		
for (ITER it = children.begin(); it != children.end(); it++)
(*it).second->push();
}
*/


 //recursively go down the tree, pushing any surplus suffixes down
void index_trie::push()
{

}



index_trie::node_type index_trie::type()
{

    if (children.empty())
    {
        switch(label->suffixes.size())
        {
        case 0: cerr<<"Shouldn't have an empty leaf."<<endl; exit(54); break;
        case 1: return uleaf; break;
        default: return leaf; break;
        }
    }
    else return branch;
}



//write a single node
void index_trie::write(cis::NUC base, 
                       ostream & treefile, 
                       ostream & posfile)
{
	
    if (is_merge_node)
    {

        std::ifstream ifile(label->treefile.c_str(), 
                            std::ios::in | std::ios::binary);
        char * itree = new char[label->treesize];
        ifile.read(itree, label->treesize);
        treefile.write(itree, label->treesize);
        if (! treefile.good())
        {
            std::fprintf(stderr, "Couldn't write %s.", 
                         label->treefile.c_str());
            exit(50);
        }
        delete itree;
        ifile.close();
        unlink(label->treefile.c_str());

        std::ifstream pfile(label->posfile.c_str(), 
                            std::ios::in | std::ios::binary);
        int64_t possize = label->poscount * POS_BYTES;
        char * ptree = new char[possize];
        pfile.read(ptree, possize);
        posfile.write(ptree, possize);
        if (! posfile.good())
        {
            std::fprintf(stderr, "Couldn't write %s.", 
                         label->posfile.c_str());
            exit(50);
        }
        delete ptree;
        pfile.close();
        unlink(label->posfile.c_str());

        return;
    }

    node_type nt = type();

    int64_t suf_bytes = longBytes(nsuffixes);
    int64_t tre_bytes = longBytes(nbytes);

    Header H;
	
    H.dat = 
    (nt<<6) 
    | (nt == branch ? tre_bytes<<3 : 0)
    | (nt == uleaf ? 0 : suf_bytes);

    assert(validHeader(H));

    switch (nt)
    {
    case leaf:
        treefile.put(H.dat);
        treefile.put(base); //is this even necessary?
        writeInt(treefile, nsuffixes, suf_bytes);
        break;

    case uleaf: {
        treefile.put(H.dat);
        treefile.put(base);
        cis::dnastring const& dsuffix = label->suffixes[0].second;

        unsigned int sbytes = dsuffix.size() / 2;
        writeInt(treefile, sbytes, 1);
        NUCpair * suffix = (NUCpair *)alloca(sbytes);
        for (unsigned int i=0; i < sbytes; i++)
        suffix[i].NN = dsuffix[2*i] | (dsuffix[2*i+1]<<4);

        treefile.write((char *)suffix, sbytes);
        break;
    }
		
    case branch: {
        treefile.put(H.dat);
        treefile.put(base);

        unsigned int ccount = 
        children.size() +
        (label->suffixes.empty() ? 0 : 1);

        writeInt(treefile, ccount, 1);
		
        writeInt(treefile, nbytes, tre_bytes);
        writeInt(treefile, nsuffixes, suf_bytes);

        //write the stub if necessary
        if (! label->suffixes.empty())
        {
            Header st;
            int nsuf = label->suffixes.size();
            int sbytes = longBytes(nsuf);
            st.dat = stub<<6 | sbytes;
            treefile.put(st.dat);
            writeInt(treefile, nsuf, sbytes);
        }
        break;
    }
		
    case stub:
        cerr<<"Shouldn't have a stub here."<<endl;
        break;
    }
	
    index_node_label::SIT sit;
	
    for (sit = label->suffixes.begin(); sit != label->suffixes.end(); sit++)
    writeInt(posfile, (*sit).first, POS_BYTES);
	
    for (ITER it = children.begin(); it != children.end(); it++)
    (*it).second->write((*it).first, treefile, posfile);

}


bool index_trie::empty()
{
    return 
    children.empty() && 
    label->suffixes.empty();
}
