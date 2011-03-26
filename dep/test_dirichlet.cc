#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <sys/timeb.h>
#include <algorithm>

int main(int argc, char ** argv)
{
    //print out the cdf of the dirichlet 
    double alpha[4];
    alpha[0] = atof(argv[1]);
    alpha[1] = atof(argv[2]);
    alpha[2] = atof(argv[3]);
    alpha[3] = atof(argv[4]);
    double cdf_cutpoint = atof(argv[5]);
    size_t num_points = static_cast<size_t>(atof(argv[6]));
    
    double theta[4];
    double cdf[4];

    double alpha0 = alpha[0] + alpha[1] + alpha[2] + alpha[3];

    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);

    for (size_t i = 0; i != 10000; ++i)
    {
        double cut = static_cast<double>(i) / 10000.0;
        printf("%10.8f\t%10.8f\n",
               cut,
               gsl_sf_beta_inc(alpha[0], alpha0 - alpha[0], cut));

    }        

    for (size_t i = 0; i != num_points; ++i)
    {
        gsl_ran_dirichlet(rand_gen, 4, alpha, theta);
        double val = 1.0;
        //double val = gsl_ran_dirichlet_pdf(4, alpha, theta);
        printf("%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n",
               theta[0], theta[1], theta[2], theta[3], val);
        
    }
    gsl_rng_free(rand_gen);

    return 0;

}
