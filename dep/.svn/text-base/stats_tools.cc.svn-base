#include "stats_tools.h"

#include <cstdio>
#include <numeric>
#include <functional>
#include <cmath>
#include <algorithm>


#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>

//calculate a new mean of n+1 samples, given the mean based on n samples
//and a new sample
double add_sample_to_mean(double mean, size_t n, double new_sample)
{
    return mean + (new_sample - mean) / static_cast<double>(n + 1);
}


double remove_sample_from_mean(double mean, size_t n, double old_sample)
{
    if (n < 2)
    {
        fprintf(stderr, "error:  remove_sample_from_mean:"
                "cannot remove a sample from a mean with only one sample.");
        exit(1);
    }

    return mean - (old_sample - mean) / static_cast<double>(n - 1);
}


double add_sample_to_covariance(double covariance, size_t n,
                                double x, double mean_x,
                                double y, double mean_y)
{
    if (n < 2)
    {
        fprintf(stderr, "error:  add_sample_to_covariance:"
                "covariance is undefined for less than two samples.");
        exit(1);
    }

    double nd = static_cast<double>(n);
    double mean_x_incr = add_sample_to_mean(mean_x, n, x);
    return ((covariance * (nd - 1))
            + ((x - mean_x_incr) * (y - mean_y)))
        / nd;
}


double remove_sample_from_covariance(double covariance, size_t n,
                                     double x, double mean_x,
                                     double y, double mean_y)
{

    if (n < 3)
    {
        fprintf(stderr, "error:  remove_sample_from_covariance:"
                "covariance is undefined for less than two samples."
                "answer will be undefined");
        exit(1);
    }

    double nd = static_cast<double>(n);
    double mean_x_decr = remove_sample_from_mean(mean_x, n, x);
    double c = covariance * (nd - 1);

    return (c - (y - mean_y) * (x - mean_x_decr)) / (nd - 2);

}



void add_to_mean_covariance_matrix(double const* point, size_t ndim, 
                                   size_t orig_num_samples,
                                   double * mean, double * covariance)
{
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            size_t d = d1 * ndim + d2;
            covariance[d] =
                add_sample_to_covariance(covariance[d], orig_num_samples,
                                         point[d1], mean[d1],
                                         point[d2], mean[d2]);
        }
    }

    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        mean[d1] = add_sample_to_mean(mean[d1], orig_num_samples, point[d1]);
    }
}



void remove_from_mean_covariance_matrix(double const* point, size_t ndim, 
                                        size_t orig_num_samples,
                                        double * mean, double * covariance)
{
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            size_t d = d1 * ndim + d2;
            covariance[d] =
                remove_sample_from_covariance(covariance[d], orig_num_samples,
                                              point[d1], mean[d1],
                                              point[d2], mean[d2]);
        }
    }

    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        mean[d1] = remove_sample_from_mean(mean[d1], orig_num_samples, 
                                           point[d1]);
    }
}


void multivariate_mean(double const* samples,
                       size_t ndim,
                       size_t num_samples,
                       double * mean)
{
    //update the mean's
    for (size_t d = 0; d != ndim; ++d)
    {
        mean[d] = gsl_stats_mean(samples + d, ndim, num_samples);
    }
}


void multivariate_mean_covariance(double const* samples,
                                  size_t ndim,
                                  size_t num_samples,
                                  double * mean, double * covariance)
{

    //update the mean's
    for (size_t d = 0; d != ndim; ++d)
    {
        mean[d] = gsl_stats_mean(samples + d, ndim, num_samples);
    }
    
    //Update the covariances
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            covariance[d1 * ndim + d2] = 
                gsl_stats_covariance_m(samples + d1, ndim,
                                       samples + d2, ndim,
                                       num_samples, mean[d1], mean[d2]);
        }
    }
}



void multivariate_autocorrelation(double const* points, 
                                  size_t ndim, 
                                  size_t num_points, size_t offset,
                                  double * autocor_matrix)
{
    
    double * mean = new double[ndim];
    multivariate_mean(points, ndim, num_points, mean);

    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        double const* points1 = points + d1;
        double sdev1 = 
            gsl_stats_sd_m(points1, ndim,
                           num_points, mean[d1]);

        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            size_t d = d1 * ndim + d2;
            double const* points2 = points + (offset * ndim) + d2;
            
            double sdev2 = gsl_stats_sd_m(points2, ndim, 
                                          num_points - offset, mean[d2]);
            double covar = 
                gsl_stats_covariance_m(points1, ndim, points2, ndim,
                                       num_points - offset, 
                                       mean[d1], mean[d2]);
                
            autocor_matrix[d] = covar / (sdev1 * sdev2);
        }
    }
    delete mean;
}



/* condense the matrix of autocorrelations to a single measure, to be
   used for overall halting criteria.
   currently, the average of diagonal elements
   !! this is clearly incorrect, but who knows what the correct measure is.
*/
double autocorrelation_measure(double const* autocor_matrix,
                               size_t ndim)
{
    double sum = 0.0;
    for (size_t d = 0; d != ndim; ++d)
    {
        sum += autocor_matrix[d * ndim + d];
    }
    return sum / static_cast<double>(ndim);
}


//returns the average measure over the offsets [1, <autocor_max_offset>]
double total_autocorrelation_measure(double const* samples,
                                     size_t ndim,
                                     size_t num_samples,
                                     size_t autocor_max_offset)
{
    double * autocor_matrix = new double[ndim * ndim];
    double * autocor_measures = new double[autocor_max_offset];

    //double mean[] = {-1, -1, -1, -1};
    for (size_t o = 1; o != autocor_max_offset; ++o)
    {
        multivariate_autocorrelation(samples, ndim,
                                     num_samples, o, autocor_matrix);

//         print_mean_covariance(stdout, mean, autocor_matrix, ndim);

        autocor_measures[o] = 
            autocorrelation_measure(autocor_matrix, ndim);
    }
    
    fprintf(stdout, "autocor:");
    for (size_t o = 1; o != autocor_max_offset; ++o)
    {
        fprintf(stdout, " %5.4g", autocor_measures[o]);
    }
    fprintf(stdout, "\n");

    double total_measure =
        std::accumulate(autocor_measures + 1, 
                        autocor_measures + autocor_max_offset,
                        0.0, std::plus<double>())
        / static_cast<double>(autocor_max_offset);
    
    if (isnan(total_measure))
    {
        total_measure = 2.0;
    }

    delete autocor_matrix;
    delete autocor_measures;
    return total_measure;
}



//returns the lowest offset between samples such that the overall autocorrelation
//at that offset is below <valid_autocor>, or if not found, returns <autocor_max_offset>
size_t best_autocorrelation_offset(double const* samples,
                                   size_t ndim,
                                   size_t num_samples,
                                   size_t autocor_max_offset,
                                   double valid_autocor)
{
    size_t msize = ndim * ndim;
    double * autocor_matrix = new double[msize];
    size_t best_offset = autocor_max_offset;
    double autocor_measure;

    //linear search to find first qualifying offset
    for (size_t o = 1; o != autocor_max_offset && o <= num_samples; ++o)
    {
        multivariate_autocorrelation(samples, ndim,
                                     num_samples, o, autocor_matrix);

        autocor_measure = autocorrelation_measure(autocor_matrix, ndim);
        if (autocor_measure < valid_autocor)
        {
            best_offset = o;
            break;
        }
    }
    delete autocor_matrix;

    return best_offset;
}


bool all_positive(double const* x, size_t n)
{
    for (size_t d = 0; d != n; ++d)
    {
        if (x[d] < 0)
        {
            return false;
        }
    }
    return true;
}

bool normalized(double const* x, size_t n, double delta)
{
    double sum = 0.0;
    for (size_t d = 0; d != n; ++d)
    {
        sum += x[d];
    }
    return gsl_fcmp(sum, 1.0, delta) == 0;
}


void normalize(double * x, size_t n, double * x_out)
{
    double sum = std::accumulate(x, x + n, 0.0);
    std::transform(x, x + n, x_out, std::bind2nd(std::divides<double>(), sum));
}




//Sample from a discrete distribution of counts, returning the partition number
int SampleDiscreteDistribution(gsl_rng * rand_gen, 
                               int const* distribution, int total_counts)
{
    int random = gsl_rng_uniform_int(rand_gen, total_counts);
    int partition;
    for (partition = 0; 1; ++partition)
    {
        if (random <= distribution[partition])
        {
            break;
        }
        random -= distribution[partition];
    }
    return partition;
}


//Sample from a discrete distribution of double values !!! check this
int SampleDiscreteDistribution(gsl_rng * rand_gen, 
                               double const* distribution, double total)
{
    //first, scale the distribution to an integer range

    double random_fraction = gsl_rng_uniform(rand_gen) * total;

    int partition;
    for (partition = 0; 1; ++partition)
    {
        if (random_fraction <= distribution[partition])
        {
            break;
        }
        random_fraction -= distribution[partition];
    }
    return partition;
}
