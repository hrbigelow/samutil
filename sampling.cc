#include "sampling.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <gsl/gsl_statistics_double.h>

#include "integrands.h"




WeightedSample::WeightedSample(size_t const _ndim, REAL const * _x, REAL const _val, REAL const _weight) : 
    ndim(_ndim), val(_val), weight(_weight)
{ 
    x = new REAL[ndim];
    cdf = new REAL[ndim];
    std::copy(_x, _x + ndim, x);
    std::fill(cdf, cdf + ndim, -1.0);
}

WeightedSample::~WeightedSample()
{
    if (x != NULL && cdf != NULL)
    {
        delete x;
        delete cdf;
    }
}

WeightedSample::WeightedSample(WeightedSample const& w)
{
    ndim = w.ndim;
    x = new REAL[w.ndim];
    std::copy(w.x, w.x + w.ndim, x);
    cdf = new REAL[w.ndim];
    std::copy(w.cdf, w.cdf + w.ndim, cdf);
    val = w.val;
    weight = w.weight;
}

WeightedSample & WeightedSample::operator=(WeightedSample const& w)
{
    if (this != &w)
    {
        ndim = w.ndim;
        x = new REAL[w.ndim];
        std::copy(w.x, w.x + w.ndim, x);
        cdf = new REAL[w.ndim];
        std::copy(w.cdf, w.cdf + w.ndim, cdf);
        val = w.val;
        weight = w.weight;
    }
    return *this;
}


// !!! check this!
bool WeightedSample::operator<(WeightedSample const& w) const
{
    if (this->ndim != w.ndim)
    {
        return this->ndim < w.ndim;
    }
    else
    {
        for (size_t d = 0; d != this->ndim; ++d)
        {
            if (this->x[d] < w.x[d])
            {
                return true;
            }
        }
        return false;
    }
}


// bool less_coordinate::operator()(COORDINATE const& c1,
//                                  COORDINATE const& c2) const
// {
    
//     if (c1.second != c2.second)
//     {
//         return c1.second < c2.second;
//     }
//     else
//     {
//         for (size_t d = 0; d != c1.second; ++d)
//         {
//             if (c1.first[d] < c2.first[d])
//             {
//                 return true;
//             }
//         }
//         return false;
//     }
// }



WeightedSample::WeightedSample() : ndim(0), x(NULL), cdf(NULL), val(0.0), weight(0.0) { }


//for sorting a 4D coordinate weighted samples vector on a given dimension
class SortCoordinate
{
    size_t dimension;
public:
    SortCoordinate(size_t dim_) : dimension(dim_) { }
    bool operator()(WeightedSample const* a,
                    WeightedSample const*  b)
    {
        return a->x[dimension] < b->x[dimension];
    }
};


class SortCDF
{
    size_t dimension;
public:
    SortCDF(size_t dim_) : dimension(dim_) { }
    bool operator()(WeightedSample const* a,
                    WeightedSample const* b)
    {
        return a->cdf[dimension] < b->cdf[dimension];
    }
};
            


class SortDimension
{
    size_t key_dimension;
public:
    SortDimension(size_t kd) : key_dimension(kd) { }
    bool operator()(double const* a, double const* b)
    {
        return a[this->key_dimension] < b[this->key_dimension];
    }
};




void find_integral_bounds(std::vector<double *> * points,
                          size_t sort_dimension,
                          double const* quantiles,
                          size_t num_quantiles,
                          double * quantile_values)
{
    std::vector<double *>::iterator start = (*points).begin();
    std::vector<double *>::iterator end = (*points).end();
    std::vector<double *>::iterator rel_start, rel_cut_point;

    size_t num_points = (*points).size();
    size_t start_offset = 0;
    size_t abs_cut_point;

    for (size_t f = 0; f != num_quantiles; ++f)
    {
        if (quantiles[f] < 0.5)
        {
            abs_cut_point = 
                std::floor(quantiles[f] * static_cast<double>(num_points));
        }
        else
        {
            abs_cut_point = 
                std::ceil(quantiles[f] * static_cast<double>(num_points));
        }
            
        rel_cut_point = start + abs_cut_point;
        rel_start = start + start_offset;

        std::nth_element(rel_start, rel_cut_point, end, SortDimension(sort_dimension));
        quantile_values[f] = rel_cut_point == end ? 0.0 : (*rel_cut_point)[sort_dimension];

        start_offset = abs_cut_point;

    }
}




void print_cdf_comparison(FILE * out_fh, 
                          AnalyticalIntegrand const* integrand,
                          std::vector<double *> * points,
                          double const* quantiles,
                          size_t const num_quantiles,
                          size_t const num_dimensions)
{

    double * quantile_buffer = new double[num_quantiles * num_dimensions];
    double ** quantile_values = new double *[num_dimensions];

    for (size_t d = 0; d != num_dimensions; ++d)
    {
        quantile_values[d] = quantile_buffer + (d * num_quantiles);
        find_integral_bounds(points, d, quantiles, 
                             num_quantiles, quantile_values[d]);
    }

    for (size_t q = 0; q != num_quantiles; ++q)
    {
//         fprintf(out_fh, "");
        for (size_t d = 0; d != num_dimensions; ++d)
        {
            
            double inverse_cdf = integrand->inv_marginal_cdf(quantiles[q], d);
            double quantile_est = integrand->marginal_cdf(quantile_values[d][q], d);
            
            fprintf(out_fh, "\t%10.8f\t%10.8f\t%10.8f\t%10.8f", 
                    quantiles[q],
                    quantile_est - quantiles[q],
                    inverse_cdf,
                    quantile_values[d][q] - inverse_cdf);
        }
        fprintf(out_fh, "\n");
        
    }
    delete quantile_buffer;
    delete quantile_values;
}


//find the mode of a 1D distribution as approximated from a set of sample points
double window_averaged_mode(std::vector<double *> * points,
                            size_t sort_dimension,
                            size_t window_size)
{
    //sort sample points
    //calculate 1/d (d being distance between point p and p + offset
    //approximate density at current point as 1/dl + 1/dr, for left and right windows.
    //at boundary (where a left or right window does not exist, substitute the density
    //at the closest point.
    size_t npoints = (*points).size();

    double * density = new double[npoints];
    std::vector<double *>::iterator start = (*points).begin();
    std::vector<double *>::iterator end = (*points).end();
    std::vector<double *>::iterator left;
    std::vector<double *>::iterator right;

    std::sort(start, end, SortDimension(sort_dimension));

    size_t index;
    
    for (left = start, right = start + window_size, index = 0; right != end; 
         ++left, ++right, ++index)
    {
        density[index] = 1.0 / ((*right)[sort_dimension] 
                                - (*left)[sort_dimension]);
    }
    //fill in the last part of density
    std::fill(density + npoints - window_size, 
              density + npoints,
              density[npoints - window_size - 1]);


    //correct density to account for left-right averaging
    for (index = npoints - 1; index != 0; --index)
    {
        size_t left_index = index > window_size ? index - window_size : 0;
        density[index] = density[index] + density[left_index];
    }

    size_t max_elem = 
        std::distance(density, 
                      std::max_element(density, density + npoints));

    delete density;
    return (*points)[max_elem][sort_dimension];

}


void print_quantiles(FILE * out_fh, 
                     std::vector<double *> * points,
                     double const* mode_point,
                     char const* line_label,
                     char const** dimension_labels,
                     char const* sums_label,
                     double const* quantiles, 
                     size_t num_quantiles,
                     size_t num_dimensions)
{

    double * quantile_sums = new double[num_quantiles];
    std::fill(quantile_sums, quantile_sums + num_quantiles, 0.0);

    double * quantile_values = new double[num_quantiles];

    double mean_sum = 0.0;
    double mode_point_sum = 0.0;

    double * mean = new double[num_dimensions];
    std::multimap<double, size_t, std::greater<double> > dim_to_mean;

    //calculate mean
    for (size_t d = 0; d != num_dimensions; ++d)
    {
        mean[d] = 0.0;
        for (size_t p = 0; p != (*points).size(); ++p)
        {
            mean[d] += (*points)[p][d];
        }
        mean[d] /= (*points).size();
        dim_to_mean.insert(std::make_pair(mean[d], d));
        mean_sum += mean[d];
    }

    //calculate mean rank order
    size_t * mean_rank_order = new size_t[num_dimensions];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm_iter;
    
    size_t d = 0;
    for (dtm_iter = dim_to_mean.begin(); dtm_iter != dim_to_mean.end(); 
         ++dtm_iter)
    {
        mean_rank_order[(*dtm_iter).second] = d;
        ++d;
    }

    for (size_t d = 0; d != num_dimensions; ++d)
    {
        fprintf(out_fh, "%s\t%s\t%Zu", line_label, dimension_labels[d], mean_rank_order[d]);
        find_integral_bounds(points, d, quantiles, num_quantiles, quantile_values);

        mode_point_sum += mode_point[d];

        fprintf(out_fh, "\t%10.8f\t%10.8f", mean[d], mode_point[d]);

        for (size_t q = 0; q != num_quantiles; ++q)
        {
            fprintf(out_fh, "\t%10.8f", quantile_values[q]);
            quantile_sums[q] += quantile_values[q];
        }
        fprintf(out_fh, "\n");
    }

    fprintf(out_fh, "%s\t%s\t%s", line_label, sums_label, sums_label);
    fprintf(out_fh, "\t%10.8f\t%10.8f", mean_sum, mode_point_sum);

    for (size_t q = 0; q != num_quantiles; ++q)
    {
        fprintf(out_fh, "\t%10.8f", quantile_sums[q]);
    }
    fprintf(out_fh, "\n");
    
    delete mean_rank_order;
    delete mean;
    delete quantile_sums;
    delete quantile_values;

}



void print_numerical_cdfs(FILE * out_fh, 
                          char const* label,
                          std::vector<double *> * points,
                          size_t ndim)
{

    std::vector<double *> sorted_coords(ndim);
    size_t stride = ndim * 2;
    double * data_buf = new double[(*points).size() * stride];
    double ** ppoints = new double *[(*points).size()];
    
    for (size_t p = 0; p != (*points).size(); ++p)
    {
        ppoints[p] = data_buf + (p * stride);
        std::copy((*points)[p], (*points)[p] + 4, ppoints[p]);
    }

    for (size_t d = 0; d != ndim; ++d)
    {
        std::sort(ppoints, ppoints + (*points).size(), SortDimension(d));
        for (size_t p = 0; p != (*points).size(); ++p)
        {
            ppoints[p][ndim+d] = p;
        }
    }

    for (size_t p = 0; p != (*points).size(); ++p)
    {
        fprintf(out_fh, "%Zu\t%s", p, label);
        for (size_t d = 0; d != ndim; ++d)
        {
            fprintf(out_fh, "\t%20.18f", ppoints[p][d]);
        }
        for (size_t d = ndim; d != stride; ++d)
        {
            fprintf(out_fh, "\t%Zu", static_cast<size_t>(ppoints[p][d]));
        }
        fprintf(out_fh, "\n");
    }
    delete data_buf;
    delete ppoints;
}



void add_normalized_dimension(double const* points, size_t ndim,
                              size_t num_points,
                              double * augmented_points_buf,
                              std::vector<double *> * augmented_points)
{
    for (size_t i = 0; i != num_points; ++i)
    {
        double const* point1 = points + (i * ndim);
        double * point2 = augmented_points_buf + (i * (ndim + 1));

        std::copy(point1, point1 + ndim, point2);
        point2[ndim] = 1.0 - std::accumulate(point2, point2 + ndim, 0.0);
        (*augmented_points)[i] = point2;
    }
}
