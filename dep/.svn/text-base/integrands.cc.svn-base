#include "integrands.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_log.h>

#include "stats_tools.h"
#include "error_estimate.h"

void Integrand::store_call(double const* x, double y)
{
    std::copy(x, x + this->ndim, this->last_point);
    this->last_value = y;
}    


bool Integrand::get_last_call(double const* x, double * y) const
{
    *y = this->last_value;
    return std::equal(this->last_point, this->last_point + this->ndim, x);
}

void ignore_underflow_handler(char const* /* reason */,
                              char const* /* file */,
                              int /* line */,
                              int /* gsl_errno */) { }



Gaussian::Gaussian(REAL const* _mean, REAL const* _covariance, 
                   size_t const _ndim, bool _is_log_form) : 
    AnalyticalIntegrand(_ndim, _is_log_form), initialized(false)
{ 
    default_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler(default_handler);

    gsl_rand_gen = gsl_rng_alloc(gsl_rng_taus);

//     mean = new REAL[ndim];
//     covariance = new REAL[ndim * ndim];
//     cdf_lowbound = new REAL[ndim];
//     cdf_highbound = new REAL[ndim];

    mean = gsl_vector_alloc(ndim);
    covariance = gsl_matrix_alloc(ndim, ndim);
    covariance_inv = gsl_matrix_alloc(ndim, ndim);
    cdf_lowbound = gsl_vector_alloc(ndim);
    cdf_highbound = gsl_vector_alloc(ndim);

    this->cholesky = gsl_matrix_alloc(this->ndim, this->ndim);

    this->set_mean(_mean);
    this->set_covariance(_covariance);
    
}


void Gaussian::Init()
{


    gsl_matrix_memcpy(this->cholesky, this->covariance);
    
    gsl_linalg_cholesky_decomp(this->cholesky);

//     gsl_matrix_set_all(this->L, 0.0);

//     for (size_t d1 = 0; d1 != ndim; ++d1)
//     {
//         for (size_t d2 = 0; d2 <= d1; ++d2)
//         {
//             double cell_value = gsl_matrix_get(this->LU, d1, d2);
//             gsl_matrix_set(this->L, d1, d2, cell_value);
//         }
//     }
    
    gsl_permutation * gp = gsl_permutation_alloc(this->ndim);
    gsl_permutation_init(gp);

    int signum;
    gsl_matrix * LU = gsl_matrix_alloc(this->ndim, this->ndim);
    gsl_matrix_memcpy(LU, this->covariance);
    gsl_linalg_LU_decomp(LU, gp, &signum);

    gsl_matrix_memcpy(this->covariance_inv, this->covariance);
    gsl_linalg_LU_invert(LU, gp, this->covariance_inv);


    double det = gsl_linalg_LU_det(LU, signum);
    gsl_matrix_free(LU);
    gsl_permutation_free(gp);

    this->factor = 1.0 / sqrt(pow(2.0 * M_PI, this->ndim) * det);
                              
//     for (size_t d = 0; d != ndim; ++d)
//     {
//         gsl_vector_set(cdf_lowbound, d, this->marginal_cdf(0.0, d));
//         gsl_vector_set(cdf_highbound, d, this->marginal_cdf(1.0, d));

//         this->normalization_constant =
//             gsl_vector_get(this->cdf_highbound, d)
//             - gsl_vector_get(this->cdf_lowbound, d);

//     }            
                              
    this->initialized = true;
}


void Gaussian::get_mean(double * _mean)
{
    for (size_t d = 0; d != this->ndim; ++d)
    {
        _mean[d] = gsl_vector_get(this->mean, d);
    }
}


void Gaussian::set_mean(double const* _mean)
{
    for (size_t d1 = 0; d1 != this->ndim; ++d1)
    {
        gsl_vector_set(mean, d1, _mean[d1]);
    }
}


void Gaussian::get_covariance(double * _covariance)
{
    for (size_t d1 = 0; d1 != this->ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != this->ndim; ++d2)
        {
            _covariance[d1 * this->ndim + d2] =
                gsl_matrix_get(this->covariance, d1, d2);
        }
    }
}


void Gaussian::set_covariance(double const* _covariance)
{
    for (size_t d1 = 0; d1 != this->ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != this->ndim; ++d2)
        {
            gsl_matrix_set(this->covariance, d1, d2, 
                           _covariance[d1 * this->ndim + d2]);
        }
    }
}


void Gaussian::update_mean(double const* _mean)
{
    this->set_mean(_mean);
    this->Init();
}


void Gaussian::update_covariance(double const* _covariance)
{
    this->set_covariance(_covariance);
    this->Init();
}


Gaussian::~Gaussian()
{
//     delete mean;
//     delete covariance;
//     delete cdf_lowbound;
//     delete cdf_highbound;
//     mean = NULL;
//     covariance = NULL;
//     cdf_lowbound = NULL;
//     cdf_highbound = NULL;

    gsl_rng_free(gsl_rand_gen);
    gsl_vector_free(mean);
    gsl_matrix_free(covariance);
    gsl_matrix_free(covariance_inv);
    gsl_matrix_free(cholesky);
    gsl_vector_free(cdf_lowbound);
    gsl_vector_free(cdf_highbound);
}


/*
  produce a sample point from a multivariate gaussian.  using cholesky
  decomposition
 */
void Gaussian::sample(double * x) const
{
    
    if (! this->initialized)
    {
        fprintf(stderr, "Gaussian not initialized.");
        exit(1);
    }

    //Piecewise sample from marginals.
    double * z = new double[this->ndim];
    for (size_t d = 0; d != this->ndim; ++d)
    {
        z[d] = gsl_ran_gaussian(this->gsl_rand_gen, 1.0);
    }

    for (size_t d1 = 0; d1 != this->ndim; ++d1)
    {
        x[d1] = gsl_vector_get(this->mean, d1);
        for (size_t d2 = 0; d2 <= d1; ++d2)
        {
            x[d1] += z[d2] * gsl_matrix_get(this->cholesky, d1, d2);
        }

    }
    delete z;
}


//gaussian in [0,1]**n, zero elsewhere
//should it be normalized?
REAL Gaussian::pdf(double const* x)
{

    if (! this->initialized)
    {
        fprintf(stderr, "Gaussian not initialized.");
        exit(1);
    }

//     static size_t recall_count = 0;
    double val;
    if (this->get_last_call(x, &val))
    {
//         recall_count++;
//         if (recall_count % 10000 == 0)
//         {
//             printf("gaussian recall count: %Zu\n", recall_count);
//         }
        return val;
    }

    gsl_set_error_handler(&ignore_underflow_handler);

    double expon;
    val = 1.0;
    
    //implement multivariate gaussian
    gsl_vector * x_offset = gsl_vector_alloc(this->ndim);
    for (size_t d = 0; d != this->ndim; ++d)
    {
        gsl_vector_set(x_offset, d, x[d]);
    }
    gsl_vector_sub(x_offset, this->mean);
    gsl_vector * partial = gsl_vector_alloc(this->ndim);
    gsl_vector_set_all(partial, 0.0);

    gsl_blas_dgemv(CblasNoTrans, 1.0, this->covariance_inv, x_offset, 0.0, partial);
    gsl_blas_ddot(x_offset, partial, &expon);

    gsl_vector_free(partial);
    gsl_vector_free(x_offset);

    val = this->factor * exp(- expon / 2.0);

    this->store_call(x, val);

    return val;
}


double Gaussian::log_pdf(double const* x)
{
    return gsl_sf_log(this->pdf(x));
}


//calculate marginal cdf according to formula
REAL Gaussian::marginal_cdf(double const xcoord, size_t const marg_dim) const
{
    double mean1 = gsl_vector_get(this->mean, marg_dim);
    double stddev = sqrt(gsl_matrix_get(this->covariance, marg_dim, marg_dim));
    double cdf = gsl_cdf_gaussian_P(xcoord - mean1, stddev);
    return cdf;
}


REAL Gaussian::inv_marginal_cdf(double const p, size_t const marg_dim) const
{

    double mean1 = gsl_vector_get(this->mean, marg_dim);
    double stddev = sqrt(gsl_matrix_get(this->covariance, marg_dim, marg_dim));

    return gsl_cdf_gaussian_Pinv(p, stddev) + mean1;
}        



void Gaussian::print_params(FILE * fh) const
{
    for (size_t d1 = 0; d1 != this->ndim; ++d1)
    {
        fprintf(fh, "%20.18g\t", gsl_vector_get(this->mean, d1));
        for (size_t d2 = 0; d2 != this->ndim; ++d2)
        {
            fprintf(fh, "\t%20.18g", gsl_matrix_get(this->covariance, d1, d2));
        }
        fprintf(fh, "\n");
    }
    
}


//currently, initializes mode point and log_scaling_factor
void Posterior::initialize(double mode_tolerance, size_t max_function_evals,
                           double const* initial_point,
                           bool verbose)
{

    this->ee->find_mode_point(1e-20, max_function_evals, 
                              initial_point,
                              this->zero_boundary,
                              verbose,
                              this->mode_point);
    
    this->log_scaling_factor = 
        this->ee->log_likelihood(this->mode_point)
        + this->ee->log_dirichlet_prior(this->mode_point);
}


REAL Posterior::pdf(double const* x)
{
    double y;
//     static size_t recall_count = 0;
    if (this->get_last_call(x, &y))
    {
        return y;
    }

    double xx[4];
    std::copy(x, x+3, xx);
    xx[3] = 1.0 - x[0] - x[1] - x[2];
    if (normalized(xx, 4, 1e-10) && all_positive(xx, 4))
    {
        y = ee->ScaledPosterior(xx, this->log_scaling_factor);
    }
    else
    {
        y = 0.0;
    }
    this->store_call(x, y);
    return y;
}


REAL Posterior::log_pdf(double const* x)
{

    double y;
//     static size_t recall_count = 0;
    if (this->get_last_call(x, &y))
    {
        return y;
    }

    double xx[4];
    std::copy(x, x+4, xx);
    if (this->ndim == 3)
    {
        xx[3] = std::max(0.0, 1.0 - x[0] - x[1] - x[2]);
    }

    if (normalized(xx, 4, 1e-10) && all_positive(xx, 4))
    {
        y = ee->log_likelihood(xx)
            + ee->log_dirichlet_prior(xx);
    }
    else
    {
        y = -DBL_MAX;
    }
    this->store_call(x, y);
    return y;
}


