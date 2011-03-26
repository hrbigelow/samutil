#include <iostream>
#include <vector>
#include <getopt.h>
#include <cstring>

#include "cisortho/dnacol.h"
#include "cisortho/dna_scanning.h"
#include "cisortho/dnastats.h"
#include "cisortho/dna.h"
#include "cisortho/string_tools.h"

std::string Outfile;
std::string DataDirectory;
std::vector<std::string> Species;
std::vector<std::string> DNA;

int usage()
{
    fprintf(stderr, 
            "Usage: make_dnas_file species_names fasta_files data_directory output_index\n\n"
            "Example: make_dnas_file hg19,rn4 hg19.cfa,rn4.cfa . genomes.dnas\n\n"
            "data_directory contains all one-liner fasta files\n"
            "output_index will contain a list of byte offsets of each contig in each fasta file\n");
    return 1;
}

int main(int argc, char ** argv)
{

    if (argc != 5)
    {
        return usage();
    }

    Species = split_token(argv[1], strlen(argv[1]), ",");
    DNA = split_token(argv[2], strlen(argv[2]), ",");
    DataDirectory = std::string(argv[3]);
    Outfile = std::string(argv[4]);

    //parse dna from SequenceFile, build a dna_collection from it
    if (Species.size() != DNA.size()) 
    {
        fprintf(stderr, "You must have the same number of species names as ");
        fprintf(stderr, "genomic sequence files.\n");
        exit(53);
    }

    cis::dna_collection dnac;
    for (size_t g=0; g < DNA.size(); g++)
    {
        try 
        {
            printf("Parsing %s genomic sequence file...\n", DNA[g].c_str());
            dnac.add(Species[g], DataDirectory, DNA[g]);
        }
        catch (string &msg)
        { 
            fprintf(stderr, "%s", msg.c_str()); return 2; 
        }
    }

    dnac.calc_offsets();
    dnac.open_dnas();

    string dnasrc(Outfile); 
    std::ofstream dnafile;
    dnafile.open(dnasrc.c_str());

    for (cis::DNAS::iterator dnac_it = dnac.begin(); dnac_it != dnac.end(); ++dnac_it) 
    dnafile<<*(*dnac_it)<<endl;
    dnafile.close();

    return 0;
}
