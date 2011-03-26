#include <iostream>
#include <fstream>
#include "index_trie.h"
#include "index_trie_scan.h"
#include "litestream.h"

//a paired down version of scan in which the node
//is merely printed out...
bool verbose = false;
//int nth_line = 2<<15;
int nth_line = 1;
bool print_bytes = true;
int num_padding_bytes = 12;
int node_header_print_width = 50;
const int MAX_PREFIX_LENGTH = 1000;

using std::cout;

void printline(node_info const& disk_node, 
               litestream & trie_file,
               string const& sprefix,
               int prefix_index,
               int counter)
{

    if (counter % nth_line == 0)
    {
        int pos_before = cout.tellp();
        if (prefix_index > 0) cout<<sprefix;
        cout<<"  "<<disk_node;
        int pos_after = cout.tellp();
        int node_padding = node_header_print_width - (pos_after - pos_before);
        for (int i=0; i < node_padding; i++) cout<<" ";

        printf("%12li%3i%12li",
               disk_node.subtree_size,
               disk_node.size,
               disk_node.suffix_count);
			

        if (print_bytes)
        {
            trie_file.seekg(- disk_node.size, litestream::POS_CURRENT);
            printf("   ");
            printBytes("-#", cout, trie_file, disk_node.size);
            int num_padding = std::max(0, ::num_padding_bytes - disk_node.size);
            for (int i=0; i < num_padding; i++)
            {
                printf("        ");
                if (i < num_padding - 1) printf(" ");
            }
        }
    }

    for (int i=0; i < disk_node.suffix_bytes; i++)
    {
        if (i == 0) cout<<"  ";
        cout<<cis::Nuc2char(disk_node.suffix[i].N1)
            <<cis::Nuc2char(disk_node.suffix[i].N2);
    }

}


void printscan(litestream & trie_file,
			   litestream & loci_file,
			   cis::dna_collection const& dnac,
			   node_info const& disk_node_parent,
			   int pi, //index of the prefix, also, zero-based depth of the trie
			   int64_t maxprint = 0,
               int maxdepth = 100)
{
	
    //bool store, parent;

    static cis::dnastring prefix(MAX_PREFIX_LENGTH, cis::M);
    static char raw_dna_chars[MAX_PREFIX_LENGTH];
    static cis::NUC raw_dna_nucs[MAX_PREFIX_LENGTH];
    static string sprefix, sraw;
    static int counter = 0;

    //cout<<cis::FromDNAString(prefix)<<'\t';
    //if ((int)prefix.size() < pi+1) prefix.resize(pi+1);
    int prefix_index_actual = disk_node_parent.type == index_trie::stub ? pi - 1 : pi;

    prefix[pi] = disk_node_parent.nuc;
    sprefix = cis::FromDNAString(prefix.substr(1, prefix_index_actual));

    //bool do_print = (counter & nth_line) == 0;
    bool do_print = true;

    if (do_print) //warning, this is changed from previous condition
    printline(disk_node_parent, trie_file, sprefix, 
              prefix_index_actual, counter);


    //only print the full line from leaves
    if (disk_node_parent.type != index_trie::branch)
    {
    
        int64_t ibytes = disk_node_parent.suffix_count * POS_BYTES;
        char * locus = new char[ibytes];
		
        loci_file.read(locus, ibytes);

        int64_t nsuffix_print = 
        maxprint != -1 ? 
        std::min(maxprint, disk_node_parent.suffix_count) : 
        disk_node_parent.suffix_count;

        int ibytes_act = nsuffix_print * POS_BYTES;

        for (int i = 0; i < ibytes_act; i+=POS_BYTES)
        {

            if (do_print)
            {
                if (i) cout<<',';
                else cout<<"  ";
            }
				
            unsigned int r = readInt(locus + i, POS_BYTES);
            std::pair<cis::dna_t const*, int64_t> dloc = dnac.FromIndex(r);
            cis::dna_t const* dna = dloc.first;
			
            dna->source().seekg(dna->seek_start_pos+dloc.second, litestream::POS_BEGIN);
            dna->source().read(raw_dna_chars, prefix_index_actual);
            cis::TranslateSeq(raw_dna_nucs, raw_dna_chars, prefix_index_actual);
            if (! strnsame(raw_dna_nucs, prefix.c_str()+1, prefix_index_actual))
            {
                cout<<cis::FromDNAString(raw_dna_nucs, prefix_index_actual)<<" != "
                    <<sprefix<<endl;

            }
			
            if (do_print)
            printf("%i", r);
            //cout<<'['<<r<<"]--"<<dloc.first->species()
            //<<":"<<dloc.first->name<<":"<<dloc.second;
        }
        delete locus;
    }
  
    int ccount = disk_node_parent.child_count;
    if (do_print) cout<<endl;
    ++counter;

    while (ccount--)
    {
        node_info disk_node_child = readNode(trie_file);
        if (pi > maxdepth)
        {
            trie_file.seekg(disk_node_child.subtree_size - disk_node_child.size,
                            litestream::POS_CURRENT);
            continue;
        }
        printscan(trie_file, loci_file, dnac, disk_node_child, pi+1, maxprint, maxdepth);
    }
}



//convert a stream of binary data into 01 representation
//usage char2bin <infile>
int main(int argc, char **argv)
{

    char const* data_directory = argv[1];
    int64_t buffer_size = 1000;
    std::ifstream trie_base;
    std::ifstream loci_base;

    litestream trie_file(trie_base, buffer_size);
    trie_file.open(argv[2]);
    litestream loci_file(loci_base, buffer_size);
    loci_file.open(argv[3]);

    std::ifstream dfile(argv[4]);

    ::print_bytes = argc >= 6 ? (strcmp(argv[5], "bytes") == 0) : false;

    int maxprint = argc == 7 ? atoi(argv[6]) : -1;

    int maxdepth = argc == 8 ? atoi(argv[7]) : 100;

    if (! data_directory)
    {
        cerr<<"Please provide a pathname to the data files."<<endl;
        exit(51);
    }

    if (! trie_file.good())
    {
        cerr<<"Couldn't open trie file "<<argv[2]<<endl;
        exit(51);
    }

    if (! loci_file.good())
    {
        cerr<<"Couldn't open gpos file "<<argv[3]<<endl;
        exit(51);
    }

    if (! dfile)
    {
        cerr<<"Couldn't open dnas file "<<argv[4]<<endl;
        exit(51);
    }

    cis::dna_collection dnac;

    dfile.seekg(0, std::ios::end);
    int endpos = dfile.tellg();
    dfile.seekg(0, std::ios::beg);

    while ((int)dfile.tellg() != endpos) 
    {
        dnac.insert(new cis::dna_t(dfile, std::string(data_directory)));
    }

    dfile.close();

    dnac.calc_offsets();
    dnac.open_dnas();

    node_info root = readNode(trie_file);

    printscan(trie_file, loci_file, dnac, root, 0, maxprint, maxdepth);

    trie_file.close();
    loci_file.close();
}
