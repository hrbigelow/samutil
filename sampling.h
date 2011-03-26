#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <utility>

#include "tools.h"

/*
  Group of functions to generate individual samples from statistical distributions

 */

class Gaussian;
class AnalyticalIntegrand;




struct WeightedSample
{
    size_t ndim;
    REAL *x;
    REAL *cdf; //uninitialized
    REAL val;
    REAL weight;

    WeightedSample(size_t const _ndim, REAL const * _x, REAL const _val, REAL const _weight);
    ~WeightedSample();
    WeightedSample(WeightedSample const& w);
    WeightedSample & operator=(WeightedSample const& w);
    WeightedSample();
    bool operator<(WeightedSample const& w) const;
    
};




typedef std::map<WeightedSample, WeightedSample> WEIGHTED_SAMPLE_MAP;


void find_integral_bounds(std::vector<double *> * points,
                          size_t sort_dimension,
                          double const* quantiles,
                          size_t num_quantiles,
                          double * quantile_values);


void print_cdf_comparison(FILE * out_fh, 
                          AnalyticalIntegrand const* integrand,
                          std::vector<double *> * points,
                          double const* quantiles,
                          size_t const num_quantiles,
                          size_t const num_dimensions);


double window_averaged_mode(std::vector<double *> * points,
                            size_t sort_dimension,
                            size_t window_size);


void print_quantiles(FILE * out_fh, 
                     std::vector<double *> * points,
                     double const* mode_point,
                     char const* line_label,
                     char const** dimension_labels,
                     char const* sums_label,
                     double const* quantiles, 
                     size_t num_quantiles,
                     size_t num_dimensions);

void print_numerical_cdfs(FILE * out_fh, 
                          char const* label,
                          std::vector<double *> * points,
                          size_t dim);


void add_normalized_dimension(double const* points, size_t ndim,
                              size_t num_points,
                              double * augmented_points_buf,
                              std::vector<double *> * augmented_points);


#endif // _SAMPLING_H
