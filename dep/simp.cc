#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <numeric>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>


#include <sys/timeb.h>

#include "sampling.h"
#include "error_estimate.h"
#include "tools.h"
#include "simulation.h"
#include "nucleotide_stats.h"
#include "base_qual_strand_reader.h"
#include "stats_tools.h"

#include "usage_strings.h"

int simp_usage()
{
    fprintf(stderr, 
            "Usage: dep simp jpd_params sim_comp > output.pileup\n"
            "sim_comp format:\n"
            "position fractionA fractionC fractionG fractionT sample_size\n\n"
);
    return 1;
}

int main_simp(int argc, char ** argv)
{

    if (argc < 3)
    {
        return simp_usage();
    }

    char * base_qual_params_file = argv[1];
    char * real_composition_file = argv[2];

    FILE * real_composition_fh = fopen(real_composition_file, "r");
    if (real_composition_fh == NULL)
    {
        fprintf(stderr, "Couldn't open real_composition file %s\n", 
                real_composition_file);
        return 0;
    }
    
    //char * founder_bases = new char[sample_size];
    
    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);

    BaseQualStrandReader reader;
    NucleotideStats model_params = reader.read_from_rdb(base_qual_params_file);

    size_t position = 0;
    size_t sample_size;

    char line[1001];

    //for each set of 4 numbers in real_base_composition, sample a founder base
    while (fgets(line, 1000, real_composition_fh) != NULL)
    {
        size_t founder_base_counts[4];
        std::fill(founder_base_counts, founder_base_counts + 4, 0.0);

        double founder_base_comp[4];

        size_t num_fields =
            sscanf(line, "%zi\t%lf\t%lf\t%lf\t%lf\t%zi\n",
                   &position,
                   founder_base_comp,
                   founder_base_comp + 1,
                   founder_base_comp + 2,
                   founder_base_comp + 3,
                   &sample_size);
        if (num_fields != 6)
        {
            fprintf(stderr, "Error: provided sample base composition file "
                    "doesn't have 6 fields per line\n");
            exit(5);
        }

        normalize(founder_base_comp, 4, founder_base_comp);

        char * called_bases = new char[sample_size];
        char * called_bases_printed = new char[sample_size];
        int * measured_quality = new int[sample_size];


        size_t majority_base_index = 
            std::distance(founder_base_comp,
                          std::max_element(founder_base_comp, founder_base_comp + 4));
        
        char fake_reference_base = Nucleotide::bases_upper[majority_base_index];

        for (size_t read_index = 0; read_index != sample_size; ++read_index)
        {
            size_t fbase_index =
                SampleDiscreteDistribution(rand_gen, founder_base_comp, 1.0);

            founder_base_counts[fbase_index]++;

            //sample a base qual strand datum from the prior distribution and print it.
            size_t datum_index =
                SampleDiscreteDistribution(rand_gen, 
                                           model_params.founder_base_likelihood[fbase_index],
                                           1.0);

            std::string const& datum_name = model_params.index_mapping[datum_index];
            BaseQualStrandReader::Datum datum = reader.get_datum_from_name(datum_name);


            //founder_bases[read_index] = Nucleotide::bases_upper[fbase_index];
            called_bases[read_index] = datum.strand == '+' ? datum.called_base : tolower(datum.called_base);
            measured_quality[read_index] = datum.quality;

            if (toupper(called_bases[read_index]) == fake_reference_base)
            {
                called_bases_printed[read_index] = datum.strand == '+' ? '.' : ',';
            }
            else
            {
                called_bases_printed[read_index] = called_bases[read_index];
            }
        }

        char infostring[100];
        sprintf(infostring,
                "founder_comp(%.5g,%.5g,%.5g,%.5g),"
                "founder_counts:(%Zu,%Zu,%Zu,%Zu)",
                founder_base_comp[0],
                founder_base_comp[1],
                founder_base_comp[2],
                founder_base_comp[3],
                founder_base_counts[0],
                founder_base_counts[1],
                founder_base_counts[2],
                founder_base_counts[3]);


        PrintPileupLine(stdout, infostring, position,
                        fake_reference_base, 
                        called_bases_printed, 
                        measured_quality, sample_size);

        delete called_bases;
        delete called_bases_printed;
        delete measured_quality;

    }
        
    //delete founder_bases;
    gsl_rng_free(rand_gen);

    fclose(real_composition_fh);
    return 0;
}
