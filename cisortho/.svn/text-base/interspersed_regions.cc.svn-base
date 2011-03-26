#include <iostream>
#include <vector>
#include <getopt.h>
#include <cstring>

#include "dna.h"
#include "dnacol.h"
#include "region.h"
#include "nested.h"

using namespace cis;

char *OutFile;
char *RegionsFile;
char *DNAIndexFile;
char * DataDirectory = const_cast<char *>(".");

int FirstRegionID;

const int64_t MAXLENGTH = 100000;

const int MaxFillLength = 100000;

using std::cout;
using std::endl;


static struct option long_options[] = {
    {"regions-file", required_argument, 0, 'i'},
    {"out-file", required_argument, 0, 'o'},
    {"dna-index-file", required_argument, 0, 'd'},
    {"data-directory", required_argument, 0, 'p'},
    {"first-region-id", required_argument, 0, 'r'},
    {0, 0, 0, 0}
};


char const* usage_message =
 "Welcome to interspersed_regions.  Retrieve all regions that are interspersed\n"
 "between the set of regions provided.  Retrieves all contiguous stretches of DNA\n"
 "(regions) NOT covered by at least one region in the set\n"
 "\n"
 "Usage: interspersed_regions [options]\n"
 "Options are as follows.  short or long options may be used\n"
 "\n"
 "-i  --regions-file     space-separated file with 'region_id' 'species' 'dna' 'start' 'end' 'strand' 'name'\n"
 "-o  --out-file         tab-separated file with 'region_id'\n"
 "-d  --dna-index-file   file containing paths and offsetse of chromosomes\n"
 "-p  --data-directory   full path to files mentioned in --dna-index-file\n"
 "-r  --first-region-id  id to assign first interpolated region.  ids count up from this value.\n";


int main(int argc, char **argv)
{

    int c;
    while (1)
    {

        int option_index = 0;
        c = getopt_long(argc, argv, "i:o:d:p:r:",
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
            
            case 'i': 
                RegionsFile = new char[strlen(optarg) + 1];
                strcpy(RegionsFile, optarg);
                break;
            case 'o': 
                OutFile = new char[strlen(optarg) + 1];
                strcpy(OutFile, optarg);
                break;
            case 'd': 
                DNAIndexFile = new char[strlen(optarg) + 1];
                strcpy(DNAIndexFile, optarg);
                break;
            case 'p':
                DataDirectory = new char[strlen(optarg) + 1];
                strcpy(DataDirectory, optarg);
                break;
            case 'r':
                FirstRegionID = atoi(optarg);
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


    cout<<"Determining DNA lengths."<<endl;

    std::string DNAIndexPath = DataDirectory;
    DNAIndexPath += std::string("/") + std::string(DNAIndexFile);

    std::ifstream species_stream(DNAIndexPath.c_str());
    if (! species_stream)
    {
        cerr<<"Couldn't open species list file "<<DNAIndexPath.c_str()<<endl;
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

    std::cout<<"Read "<<dnac.size()<<" pieces of DNA, "
             <<dnac.num_bases()<<" total bases."<<endl;

  
    cout<<"Retrieving regions from RDB source "<<RegionsFile<<endl;

    FILE* out_stream = fopen(OutFile, "a");
    if (NULL == out_stream)
    {
        std::fprintf(stderr, "Couldn't open output sequence file %s\n", OutFile);
        exit(56);
    }

    region::region_id = FirstRegionID;


    bool without_group = false;
    cis::REG_MAP RegionsToRetrieve = RDBToRegions(dnac, RegionsFile, without_group);

    for (cis::REG_MAP::iterator p = RegionsToRetrieve.begin(); 
         p != RegionsToRetrieve.end(); ++p)
    {

        cis::dna_t const* dna = (*p).first;
        REGIONS_MULTI & regions = (*p).second;

        REG_PC start_region =
        new region(*dna, std::string(""), 0, 0, POS, -1, -1, region::NONE);
        regions.insert(start_region);

        REG_PC end_region =
        new region(*dna, std::string(""), dna->length(), dna->length(), POS, 
                   -1, -1, region::NONE);

        regions.insert(end_region);

        region_tree * rtree = BuildTree(regions);
     
        REG_PC prev_reg = NULL;
        REG_PC curr_reg = NULL;

        REGIONS_MULTI intervening_regions;
        REG_PC intervening_region;

        cis::BY_START::iterator node_it;

        for (node_it = rtree->children->begin(); 
             node_it != rtree->children->end(); ++node_it)
        {
            prev_reg = curr_reg;
            curr_reg = (*node_it).second->interval;

            if (prev_reg == NULL || 
                RegionDistance(*prev_reg, *curr_reg) < 1) continue;

            //create intervening region on positive strand
            intervening_region = 
            new region(*dna, std::string(""), 
                       prev_reg->end, curr_reg->start,
                       POS, -1, -1, region::NONE);

            intervening_regions.insert(intervening_region);

            //create intervening region on negative strand
            intervening_region = 
            new region(*dna, std::string(""),
                       prev_reg->end, curr_reg->start,
                       NEG, -1, -1, region::NONE);

            intervening_regions.insert(intervening_region);

        }
   

        for (RIT rit = intervening_regions.begin(); 
             rit != intervening_regions.end(); ++rit)
        {
            REG_PC dna_region = (*rit);
            PrintDNARegion(dna_region, out_stream);
            delete dna_region;

        }

        delete rtree;
   
    }

    fclose(out_stream);
    return 0;
}
