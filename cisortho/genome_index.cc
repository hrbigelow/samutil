#include <iostream>
#include <sstream>
#include <cstring>

#include <getopt.h>

#include "index_trie.h"
#include "dnacol.h"
#include "dna_scanning.h"
#include "dnastats.h"
#include "dna.h"
#include "litestream.h"

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>


using std::cerr;
using std::cout;
using std::endl;
using std::string;

int MinTreeDepth = 6;
int MaxTreeDepth = 50;
int MaxSuffixLoad = 2;
int PrefixLength = 2;
int ChunkSize = 100000;
char *OutfilePrefix = NULL;
char *SinglePrefix = NULL;
char *FolderPrefix = NULL;
char *SpeciesNames = NULL;
char *SequenceFiles = NULL;

char * DataDirectory = const_cast<char *>(".");

int n_DNA, n_Species;
char **DNA, **Species;


static struct option long_options[] = {
    {"min-tree-depth", required_argument, 0, 'x'},
    {"max-tree-depth", required_argument, 0, 'y'},
    {"prefix-length", required_argument, 0, 'p'},
    {"max-suffix-load", required_argument, 0, 'l'},
    {"chunk-size", required_argument, 0, 'c'},
    {"species-names", required_argument, 0, 'n'},
    {"dna-directory", required_argument, 0, 't'},
    {"outfile-prefix", required_argument, 0, 'o'},
    {"single-prefix", required_argument, 0, 'f'},
    {"folder-prefix", required_argument, 0, 'd'},
    {0, 0, 0, 0}
};


char const* usage_message =
 "Usage: genome_index [options]\n"
 "Options are as follows.  short or long options may be used\n"
 "\n"
 "-x  --min-tree-depth  minimum pattern length to search.\n"
 "-y  --max-tree-depth  maximum depth of the index tree.\n"
 "-p  --prefix-length  length of prefix (3 is good) for splitting the subtree building.\n"
 "-l  --max-suffix-load  maximum number of suffixes at a leaf before pushing further\n"
 "-c  --chunk-size  size of reading a chunk of DNA sequence for scanning\n"
 "-n  --species-names  species names (comma-separated)\n"
 "-s  --sequence-files  sequence files (comma-separated)\n"
 "-t  --dna-directory  directory containing 'sequence-files'\n"
 "-o  --outfile-prefix  prefix on output files with extensions 'gpos', 'itree', and 'dnas'\n"
 "-f  --single-prefix  for testing\n"
 "-d  --folder-prefix  for partial runs\n";


int main(int argc, char ** argv)
{
    int c;
    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "x:y:p:l:c:n:t:o:f:d:",
                        long_options, &option_index);
    
        if (c == -1) 
        {
            printf("%s\n", usage_message);
            return 0;
            break;

            switch (c)
            {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                printf (" with arg %s", optarg);
                printf ("\n");
                break;
            
            case 'x': 
                MinTreeDepth = atoi(optarg); 
                break;
            case 'y': 
                MaxTreeDepth = atoi(optarg);
                break;
            case 'p': 
                PrefixLength = atoi(optarg);
                break;
            case 'l':
                MaxSuffixLoad = atoi(optarg);
                break;
            case 'c':
                ChunkSize = atoi(optarg);
                break;
            case 'n': 
                SpeciesNames = new char[strlen(optarg) + 1];
                strcpy(SpeciesNames, optarg);
                break;
            case 's': 
                SequenceFiles = new char[strlen(optarg) + 1];
                strcpy(SequenceFiles, optarg);
                break;
            case 't': 
                DataDirectory = new char[strlen(optarg) + 1];
                strcpy(DataDirectory, optarg);
                break;
            case 'o':
                OutfilePrefix = new char[strlen(optarg) + 1];
                strcpy(OutfilePrefix, optarg);
                break;
            case 'f':
                SinglePrefix = new char[strlen(optarg) + 1];
                strcpy(SinglePrefix, optarg);
                break;
            case 'd':
                FolderPrefix = new char[strlen(optarg) + 1];
                strcpy(OutfilePrefix, optarg);
                break;


                printf ("option -f with value `%s'\n", optarg);
                break;
     
            case '?':
                /* getopt_long already printed an error message. */
                break;
     
            default:
                abort ();
            }
        }
    }

    if (OutfilePrefix == NULL)
    {
        cerr<<"Please provide an output file prefix (option -o)"<<endl;
        exit(53);
    }

	
    std::ostringstream TreeTempDir; 
    std::ostringstream PosTempDir; 

    if (FolderPrefix == NULL)
    {
        PosTempDir<<getpid()<<"_gpos";
        TreeTempDir<<getpid()<<"_itrie";
    }
    else 
    {
        PosTempDir<<FolderPrefix<<"_gpos";
        TreeTempDir<<FolderPrefix<<"_itrie";
    }
	
    //parse dna from SequenceFile, build a dna_collection from it
    if (n_Species != n_DNA) 
    {
        cerr<<"You must have the same number of species names as genomic sequence files."<<endl;
        exit(53);
    }

    std::vector<string> sSpecies(n_Species);
    std::string sDataDirectory(DataDirectory);

    for (int g=0; g < n_Species; g++) sSpecies[g] = Species[g];

    cis::dna_collection dnac;
    for (int g=0; g < n_DNA; g++)
    {
        try 
        {
            std::cout<<"Parsing "<<DNA[g]<<" genomic sequence file..."<<std::endl;
            dnac.add(sSpecies[g], sDataDirectory, DNA[g]);
        }
        catch (string &msg)
        { 
            cerr<<msg<<endl; return 2; 
        }
    }

    dnac.calc_offsets();
    dnac.open_dnas();

    mkdir(TreeTempDir.str().c_str(), 0755);
    mkdir(PosTempDir.str().c_str(), 0755);
    DIR * treedir = opendir(TreeTempDir.str().c_str());
    DIR * posdir = opendir(PosTempDir.str().c_str());

    std::ofstream treefile, posfile;
    std::ifstream tree_bufbase;
    litestream tree_bufstream(tree_bufbase, ChunkSize);

    index_trie::set_max_suffix_load(MaxSuffixLoad);
    index_trie * root = new index_trie();
  
    cout<<"Scanning all dna for prefixes...";
    cout.flush();

    std::set<cis::dnastring> Prefixes, ptmp;
    std::set<cis::dnastring>::iterator prefix_iter;

    if (FolderPrefix != NULL)
    {
        //partial run.  build a tree from already run subtrees in a subfolder
        dirent * current_file;
        while ((current_file = readdir(treedir)) != NULL)
        {

            string prefix = current_file->d_name;
            //string prefix = itr->leaf();
            cis::dnastring dprefix = cis::ToDNAString(prefix);
            Prefixes.insert(dprefix);
            
            cis::dnastring pfx;
            for (prefix_iter = Prefixes.begin();
                 prefix_iter != Prefixes.end(); ++prefix_iter)
            {
                
                pfx = *prefix_iter;
                
                string spfx = cis::FromDNAString(pfx);
                
                string treepath = TreeTempDir.str() + spfx;
                string pospath = PosTempDir.str() + spfx;
                
                tree_bufstream.open(treepath.c_str());
                
                if (! tree_bufstream.good())
                {
                    cerr<<"Couldn't open "<<spfx<<" in "<<treepath<<endl;
                    exit(50);
                }

                node_info subroot = readNode(tree_bufstream);

                root->insertfile(pfx.c_str(), pfx.size(),
                                 treepath, subroot.subtree_size,
                                 pospath, subroot.suffix_count);

                tree_bufstream.close();
            }
        }
    }

    else 
    {
        for (cis::DNAS::iterator dna_it = dnac.begin(); 
             dna_it != dnac.end(); ++dna_it)
        {
      
            ptmp = dnastats::Spectrum(*(*dna_it), PrefixLength, ChunkSize);
            Prefixes.insert(ptmp.begin(), ptmp.end());
        }
    
        cout<<"found "<<Prefixes.size()<<" prefixes."<<endl;
	
        cis::dnastring pfx;
        for (prefix_iter = Prefixes.begin();
             prefix_iter != Prefixes.end(); ++prefix_iter)
        {
      
            pfx = *prefix_iter;
      
            index_trie * subroot = new index_trie();
            subroot->build_subtree(pfx, dnac, MinTreeDepth, MaxTreeDepth, ChunkSize);
            subroot->push();
      
            cis::NUC lastNuc = *pfx.rbegin();
      
            if (! subroot->empty())
            {
        
                subroot->count();
        
                string spfx = cis::FromDNAString(pfx);
                string treepath = TreeTempDir.str() + string("/") + spfx;
                string pospath = PosTempDir.str() + string("/") + spfx;
		
                treefile.open(treepath.c_str(), std::ios::out | std::ios::binary);
                if (! treefile)
                {
                    cerr<<"Couldn't open "<<spfx<<" in "<<treepath<<endl;
                    exit(50);
                }
        
                posfile.open(pospath.c_str(), std::ios::out | std::ios::binary);
                if (! posfile)
                {
                    cerr<<"Couldn't open "<<spfx<<" in "<<treepath<<endl;
                    exit(50);
                }

                cout<<"Writing subtree with prefix "<<spfx<<", "
                    <<subroot->nbytes<<" bytes, "
                    <<subroot->nsuffixes<<" patterns...";
                cout.flush();
		
                subroot->write(lastNuc, treefile, posfile);
                treefile.close();
                posfile.close();

                root->insertfile(pfx.c_str(), pfx.size(), 
                                 treepath, subroot->nbytes,
                                 pospath, subroot->nsuffixes);

                cout<<"done."<<endl;

                //return 0;

            } 
		
            delete subroot; //memory
        }

    }



    string tree(OutfilePrefix); tree += ".itrie";
    treefile.open(tree.c_str(), std::ios::out | std::ios::binary);

    string pos(OutfilePrefix); pos += ".gpos";
    posfile.open(pos.c_str(), std::ios::out | std::ios::binary);

    //return 0;

    root->push();
    root->count();
    root->write(cis::X, treefile, posfile);

    cout<<"Wrote "<<root->nbytes<<" bytes and "
        <<root->nsuffixes<<" patterns."<<endl;

    delete root;

    treefile.close();
    posfile.close();

    string dnasrc(OutfilePrefix); dnasrc += ".dnas";
    std::ofstream dnafile;
    dnafile.open(dnasrc.c_str());

    for (cis::DNAS::iterator dna_it = dnac.begin(); 
         dna_it != dnac.end(); ++dna_it) 
    dnafile<<*(*dna_it)<<endl;

    dnafile.close();


    //delete all files in treedir, posdir
    cout<<"Removing all files in "
        <<TreeTempDir.str().c_str()
        <<" and "<<PosTempDir.str().c_str()<<"."<<endl;

    dirent * current_file;
    while ((current_file = readdir(treedir)) != NULL) 
    {
        cout<<"Would be unlinking "<<current_file->d_name<<endl;
        //unlink(current_file->d_name);
    }
    rmdir(TreeTempDir.str().c_str());

    while ((current_file = readdir(posdir)) != NULL) 
    {
        cout<<"Would be unlinking "<<current_file->d_name<<endl;
        //unlink(current_file->d_name);
    }
    rmdir(PosTempDir.str().c_str());
  
    return 0;
}
