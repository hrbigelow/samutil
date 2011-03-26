#include "nucleotide_reader.h"
#include "nucleotide_stats.h"

#include <cstdlib>

NucleotideReader::NucleotideReader() : fh(NULL)
{
}

NucleotideReader::~NucleotideReader()
{
    if (fh != NULL)
    {
        fclose(fh);
    }
}


bool NucleotideReader::initialize(char const* file)
{
    fh = fopen(file, "r");
    if (fh == NULL)
    {
        fprintf(stderr, "NucleotideReader: couldn't open nucleotide file %s\n",
                file);
        return false;
    }
    return true;
}


bool NucleotideReader::more_loci()
{
    bool more = true;
    int next_char = fgetc(this->fh);
    if (next_char == EOF)
    {
        more = false;
    }
    ungetc(next_char, this->fh);
    return more;
}



JPD_DATA parse_jpd_rdb_file(char const* rdb_file);


NucleotideStats NucleotideReader::read_from_rdb(char const* rdb_file)
{

    JPD_DATA raw_counts = parse_jpd_rdb_file(rdb_file);
    NucleotideStats stats(raw_counts.size());
    stats.initialize(raw_counts);

    return stats;
}



JPD_DATA parse_jpd_rdb_file(char const* rdb_file)
{
    //intialize data_prior
    FILE * rdb_fh = fopen(rdb_file, "r");
    
    if (rdb_fh == NULL)
    {
        fprintf(stderr, "Couldn't open data jpd file %s\n", 
                rdb_file);
        exit(5);
    }

    JPD_DATA counts_map;
    char name[100];
    double counts[4];
    double counts_sum;

    while (! feof(rdb_fh))
    {
        fscanf(rdb_fh, "%s\t%lf\t%lf\t%lf\t%lf\n", name, counts, counts+1, counts+2, counts+3);

        for (size_t bi = 0; bi != 4; ++bi)
        {
            if (counts[bi] < 0)
            {
                fprintf(stderr, "BaseQualStrandReader::parse_rdb_file: "
                        "found negative count for datum %s.\n", name);
                exit(11);
            }
        }
        counts_sum = counts[0] + counts[1] + counts[2] + counts[3];

        if (counts_sum == 0)
        {
            continue;
        }

        if (counts_map.find(std::string(name)) != counts_map.end())
        {
            fprintf(stderr, "Encountered duplicate data name: %s", name);
            exit(8);
        }
        counts_map.insert(std::make_pair(std::string(name), nuc_frequency(counts)));
    }
    fclose(rdb_fh);

    return counts_map;
}
