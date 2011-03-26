#ifndef _STATS_TOOLS_H
#define _STATS_TOOLS_H

#include <cstdlib>
#include <gsl/gsl_rng.h>


/*
  calculate <ndim>-dimensional mean from the samples packed as
  [x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]
  but only of every <sample_stride> -th sample.
  <num_samples> is the total number of samples, regardless of how many
  are used to calculate the mean
*/
void multivariate_mean(double const* samples,
                       size_t ndim,
                       size_t num_samples,
                       double * mean);


//calculate the mean and covariance from the set of samples packed as
//[x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]
void multivariate_mean_covariance(double const* samples,
                                  size_t ndim,
                                  size_t num_samples,
                                  double * mean, double * covariance);


//calculate the mean and covariance from the set of samples packed as
//[x1, y1, z1, x2, y2, z2, x3, y3, z3, ... ]
void multivariate_mean_covariance(double const* samples,
                                  size_t ndim,
                                  size_t num_samples,
                                  double * mean, double * covariance);


/*
  calculate multivariate autocorrelation on the sequence of <points>,
  assumed to be of dimension <ndim> and <num_points> at the offset
  given.  take coordinates mod <num_points> (wrap around)
 */
void multivariate_autocorrelation(double const* points, size_t ndim, 
                                  size_t num_points, size_t offset,
                                  double * autocorrelation);


double add_sample_to_mean(double mean, size_t n, double new_sample);

double remove_sample_from_mean(double mean, size_t n, double old_sample);

/*
  calculate covariance of n+1 samples.  given: covariances of n
  samples, a new sample (x, y) to add, and the mean (mean_x, mean_y)
  of the n samples before adding the new one.
*/
double add_sample_to_covariance(double covariance, size_t n,
                                double x, double mean_x,
                                double y, double mean_y);


/* 
   calculate the covariance of n-1 samples.  given: covariance of n
   samples, the pair of samples (x, y) to be removed, and the means
   (mean_x, mean_y) of the n samples before removal.
*/
double remove_sample_from_covariance(double covariance, size_t n,
                                     double x, double mean_x,
                                     double y, double mean_y);



/*
  update an <ndim>-dimensional <mean> vector and <covariance> matrix,
  originally calculated from <orig_num_samples> samples, with a new
  <point>.
 */
void add_to_mean_covariance_matrix(double const* point, size_t ndim, 
                                   size_t orig_num_samples,
                                   double * mean, double * covariance);



/*
  update an <ndim>-dimensional <mean> vector and <covariance> matrix,
  originally calculated from <orig_num_samples>, by removing the given
  <point>, assumed to be one of the original sample points used for
  the mean and covariance as provided.
 */
void remove_from_mean_covariance_matrix(double const* point, size_t ndim, 
                                        size_t orig_num_samples,
                                        double * mean, double * covariance);



double autocorrelation_measure(double const* autocor_matrix,
                               size_t ndim);


double total_autocorrelation_measure(double const* samples,
                                     size_t ndim,
                                     size_t num_samples,
                                     size_t autocor_max_offset);


size_t best_autocorrelation_offset(double const* samples,
                                   size_t ndim,
                                   size_t num_samples,
                                   size_t autocor_max_offset,
                                   double valid_autocor);

bool all_positive(double const* x, size_t n);

bool normalized(double const* x, size_t n, double delta);


void normalize(double * x, size_t n, double * x_out);


//Sample from a discrete distribution of counts, returning the partition number
int SampleDiscreteDistribution(gsl_rng * rand_gen,
                               int const* distribution, int total_counts);

//Sample from a discrete distribution of double values
int SampleDiscreteDistribution(gsl_rng * rand_gen,
                               double const* distribution, double total);


#endif // _STATS_TOOLS_H
