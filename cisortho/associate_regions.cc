#include <cmath>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <getopt.h>
#include <cstring>

#include "string_tools.h"
#include "cluster.h"
#include "nested.h"
#include "region_association.h"
#include "dnacol.h"

//Read PWM file into a vector<*int>
//(since our bit accuracy is limited, we don't need double
//and PWMsize

//using namespace cis;

using std::cout;

char *DNAIndexFile = NULL;
char *DataDirectory = NULL;

char *TargetRegionFile;
char *QueryRegionFile;
char *OutFile;

int MaxAssociationDist = 10000;
int MinAssociationDist = -10000; 
int MaxRegionsBetween = 2;

static struct option long_options[] = {
  {"max-association-dist", required_argument, 0, 'b'},
  {"min-association-dist", required_argument, 0, 'c'},
  {"max-regions-between", required_argument, 0, 'm'},
  {"out-file", required_argument, 0, 'o'},
  {"dna-index-file", required_argument, 0, 'd'},
  {"data-directory", required_argument, 0, 'p'},
  {"query-regions-file", required_argument, 0, 'q'},
  {"target-region-file", required_argument, 0, 't'},
  {0, 0, 0, 0}
};


int main (int argc, char ** argv) {

  int c;

  while (1) {

    int option_index = 0;
    c = getopt_long(argc, argv, "b:c:m:x:o:d:p:q:t:",
                    long_options, &option_index);
    
    if (c == -1) break;

    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf ("option %s", long_options[option_index].name);
      if (optarg)
        printf (" with arg %s", optarg);
      printf ("\n");
      break;
      
    case 'b': MaxAssociationDist = atoi(optarg); break;
    case 'c': MinAssociationDist = atoi(optarg); break;
    case 'm': MaxRegionsBetween = atoi(optarg); break;
    case 'o': OutFile = new char[strlen(optarg) + 1]; 
      strcpy(OutFile, optarg); break;
    case 'd': DNAIndexFile = new char[strlen(optarg) + 1];
      strcpy(DNAIndexFile, optarg); break;
    case 'p': DataDirectory = new char[strlen(optarg) + 1];
      strcpy(DataDirectory, optarg); break;
    case 'q': QueryRegionFile = new char[strlen(optarg) + 1];
      strcpy(QueryRegionFile, optarg); break;
    case 't': TargetRegionFile = new char[strlen(optarg) + 1];
      strcpy(TargetRegionFile, optarg); break;


      printf ("option -f with value `%s'\n", optarg);
      break;
     
    case '?':
      /* getopt_long already printed an error message. */
      break;
     
    default:
      abort ();
    }
  }
     
    
  /* Print any remaining command line arguments (not options). */
  if (optind < argc) {
      printf("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
  }
     
   FILE * outstream = fopen(OutFile, "w");
  if (outstream == NULL){
    fprintf(stderr,
            "Can't open region relations output file '%s'.",
            OutFile);
    exit(50);
  }
  

  char OrigDir[1000];
  getcwd(OrigDir, 1000);

  cout<<"Determining DNA lengths."<<endl;
	
  //   std::string DNAIndexPath(DataDirectory);
  //   DNAIndexPath += std::string("/") + std::string(DNAIndexFile);

  std::ifstream DNAIndex_stream(DNAIndexFile);
  if (! DNAIndex_stream){
    cerr<<"Couldn't open species list file "<<DNAIndexFile<<endl;
    exit(50);
  }

  cis::dna_collection dnac;

  DNAIndex_stream.seekg(0, std::ios::end);
  int endpos = DNAIndex_stream.tellg();
  DNAIndex_stream.seekg(0, std::ios::beg);

  while ((int)DNAIndex_stream.tellg() != endpos) 
    dnac.insert(new cis::dna_t(DNAIndex_stream, std::string(DataDirectory)));
  DNAIndex_stream.close();
  dnac.calc_offsets();


  bool without_group = false;
  cout<<"Loading target regions from "<<TargetRegionFile<<"...";
  cout.flush();
  cis::REG_MAP target_regions_dna_map = 
    RDBToRegions(dnac, TargetRegionFile, without_group);
  cout<<"Done."<<endl;


  //now parse HSP regions, which are a generalization of orthologs...
  cout<<"Loading regions to associate, from "<<QueryRegionFile<<"...";
  cout.flush();
  cis::REG_MAP query_regions_dna_map = 
    RDBToRegions(dnac, QueryRegionFile, without_group);
  cout<<"Done."<<endl;

  int last_id = PrintRegionAssociations(query_regions_dna_map, target_regions_dna_map,
                                        MaxAssociationDist, MinAssociationDist, 
                                        MaxRegionsBetween, 0, outstream);

  printf("Wrote %i region associations.\n", last_id);

  fclose(outstream);

  delete OutFile;
  delete DNAIndexFile;
  for (cis::DNAS::iterator dnacol_iter = dnac.begin();
       dnacol_iter != dnac.end(); ++dnacol_iter)
    delete (*dnacol_iter);
  
  return 0;

}
