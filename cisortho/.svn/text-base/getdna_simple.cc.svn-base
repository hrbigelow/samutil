#include <iostream>
#include <vector>
#include <getopt.h>
#include <cstring>

#include "dna.h"
#include "dnacol.h"
#include "region.h"
//#include "nested.h"

using namespace cis;

char *OutFile;
char *OutDNAFile;

char *QueryRegionsFile;
char *DNAIndexFile;
char *DataDirectory;

const int64_t MAXLENGTH = 100000;

const int MaxFillLength = 100000;

using std::cout;
using std::endl;

static struct option long_options[] = {
    {"query-regions-file", required_argument, 0, 'i'},
    {"out-region-file", required_argument, 0, 'o'},
    {"out-dna-file", required_argument, 0, 's'},
    {"dna-index-file", required_argument, 0, 'd'},
    {"data-directory", required_argument, 0, 'p'},
    {0,0,0,0}
};


int main(int argc, char ** argv)
{

    int c;

    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "i:o:s:d:p:", long_options, &option_index);
    
        if (c == -1) break;

        switch (c)
        {
        case 0:
            if (long_options[option_index].flag != 0) break;
            printf ("option %s", long_options[option_index].name);
            if (optarg) printf (" with arg %s\n", optarg);
            break;
      
        case 'i': QueryRegionsFile = new char[strlen(optarg) +1]; 
            strcpy(QueryRegionsFile, optarg); break;

        case 'o': OutFile = new char[strlen(optarg) + 1]; 
            strcpy(OutFile, optarg); break;

        case 's': OutDNAFile = new char[strlen(optarg) + 1]; 
            strcpy(OutDNAFile, optarg); break;

        case 'd': DNAIndexFile = new char[strlen(optarg) + 1]; 
            strcpy(DNAIndexFile, optarg); break;

        case 'p': DataDirectory = new char[strlen(optarg) + 1]; 
            strcpy(DataDirectory, optarg); break;

        default:
            printf("usage: getdna_simple -i <query_region_file> -o <out_region_file>\n"
                   "-s <out_dna_file> -d <dna_index_file> -p <data_directory>\n"
                   "\n"
                   "<query_region_file>:  rdb file with regions to retrieve with fields\n"
                   "(region_id group_id species dna start end strand name)\n"
                   "<out_region_file>: output file with (region_id, species, dna, start, end, strand)\n"
                   "<out_dna_file>: output file with (region_id, sequence)\n"
                   "<dna_index_file>: file containing paths and offsets of chromosomes\n"
                   "<data_directory>: full path to files mentioned in species-file"
                   );
            exit(1);
            break;
        }
    }

    printf("Determining DNA lengths.");
  
    std::ifstream species_stream(DNAIndexFile);
    if (! species_stream)
    {
        cerr<<"Couldn't open species list file "<<DNAIndexFile<<endl;
        exit(50);
    }

    cis::dna_collection dnac;

    species_stream.seekg(0, std::ios::end);
    int endpos = species_stream.tellg();
    species_stream.seekg(0, std::ios::beg);

    while ((int)species_stream.tellg() != endpos) 
    dnac.insert(new cis::dna_t(species_stream, DataDirectory));
    species_stream.close();

    dnac.calc_offsets();
    dnac.open_dnas();

    std::cout<<"Read "<<dnac.size()<<" pieces of DNA, "<<dnac.num_bases()<<" total bases."<<endl;

    //construct genes
    //these are assumed exon annotations
    cout<<"Parsing regions to retrieve dna from tab-separated file: "<<endl
        <<QueryRegionsFile<<endl;

    bool without_group = false;
    cis::REG_MAP ToRetrieve = RDBToRegions(dnac, QueryRegionsFile, without_group);


    FILE* outstream1 = fopen(OutFile, "w");
    if (NULL == outstream1)
    {
        std::fprintf(stderr, "Couldn't open output sequence file %s\n", OutFile);
        exit(56);
    }

    FILE* outstream2 = fopen(OutDNAFile, "w");
    if (NULL == outstream2)
    {
        std::fprintf(stderr, "Couldn't open output sequence file %s\n", OutFile);
        exit(56);
    }
  

    for (cis::REG_MAP::iterator p = ToRetrieve.begin(); p != ToRetrieve.end(); ++p)
    {
        REGIONS_MULTI & regions = (*p).second;
     
        for (RIT rit = regions.begin(); rit != regions.end(); ++rit)
        {
            PrintDNARegion(*rit, outstream1);
            PrintDNARegionSimple(*rit, ::MaxFillLength, outstream2);
        }

        for (RIT rit = regions.begin(); rit != regions.end(); ++rit) 
        delete (*rit);
    }
    fclose(outstream1);
    fclose(outstream2);

    return 0;
}
