#include <cstdio>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <sys/timeb.h>

#include "stats_tools.h"
#include "simulation.h"


int simc_usage()
{
    fprintf(stderr, 
            "Usage: dep simc [options] basecomp_priors nsim > sim_comp\n"
            "Options:\n\n"
            "-d INT     sequencing depth\n\n"
            "\n"
            "basecomp_priors format:\n"
            "fractionA fractionC fractionG fractionT locus_fraction_or_count\n\n"
            "output format:\n"
            "locus_number fractionA fractionC fractionG fractionT sequence_depth\n\n"
            );
    return 1;
}

int main_simc(int argc, char ** argv)
{
    if (argc < 3)
    {
        return simc_usage();
    }

    size_t sequence_depth = 100;
    char c;

    while ((c = getopt(argc, argv, "d:")) >= 0)
    {
        switch(c)
        {
        case 'd': sequence_depth = static_cast<size_t>(atof(optarg)); break;
        default: return simc_usage(); break;
        }
    }
    if (argc - optind < 2)
    {
        return simc_usage();
    }

    char * basecomp_prior_file = argv[optind];
    size_t number_lines = static_cast<size_t>(atof(argv[optind + 1]));

    //read the whole file

    
    double * basecomp;
    double * basecomp_dist;

    double * numbers_buffer = 
        parse_basecomp_prior_file(basecomp_prior_file,
                                  &basecomp, &basecomp_dist);

    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);

    for (size_t l = 0; l != number_lines; ++l)
    {
        size_t basecomp_index = 
            SampleDiscreteDistribution(rand_gen, basecomp_dist, 1.0);

        double * p = basecomp + (4 * basecomp_index);
        fprintf(stdout, "%Zu\t%12.10f\t%12.10f\t%12.10f\t%12.10f\t%Zu\n", 
                l, p[0], p[1], p[2], p[3], sequence_depth);
    }
                
    gsl_rng_free(rand_gen);
    delete numbers_buffer;

    return 0;
}
