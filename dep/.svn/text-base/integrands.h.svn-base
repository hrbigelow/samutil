#ifndef _INTEGRANDS_H
#define _INTEGRANDS_H

#include "sampling.h"
#include "tools.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>


class Integrand
{
protected:
    size_t ndim;
    
private:
    double * last_point;
    double last_value;
    
    
public:

    bool may_underflow;

    Integrand(size_t _ndim, bool _may_underflow) : 
        ndim(_ndim), may_underflow(_may_underflow)
    { 
        last_point = new double[ndim];
    }
    ~Integrand()
    {
        delete last_point;
    }
    void store_call(double const* x, double y);
    bool get_last_call(double const* x, double * y) const;

    virtual double pdf(double const* x) = 0;
    virtual double log_pdf(double const* x) = 0;

};


class SamplingFunction : public Integrand
{
public:
    SamplingFunction(size_t _ndim, bool _may_underflow) : Integrand(_ndim, _may_underflow) { }
    virtual void sample(double * x) const = 0;
    virtual void sample_conditioned(double const* x_tau,
                                    double * x_star) = 0;
};


/*
  want a class that allows to sample, evaluate, memoize.
 */


class AnalyticalIntegrand : public SamplingFunction
{
public:
    AnalyticalIntegrand(size_t ndim, bool _may_underflow) : 
        SamplingFunction(ndim, _may_underflow) { }

    virtual double marginal_cdf(double const xcoord, size_t const marg_dim) const = 0;
    virtual double inv_marginal_cdf(double const p, size_t const marg_dim) const = 0;
};

/*
  implementation of multivariate gaussian
 */
class Gaussian : public AnalyticalIntegrand
{

    gsl_vector * mean, * cdf_lowbound, * cdf_highbound;
    gsl_matrix * covariance, * covariance_inv, * cholesky, * L;

/*     size_t ndim; */
    
    double factor;
    double normalization_constant;

    gsl_rng * gsl_rand_gen;
    gsl_error_handler_t * default_handler;

    bool initialized;

public:
    
    Gaussian(REAL const* _mean, REAL const* _covariance, size_t const _ndim,
             bool _may_underflow);

    ~Gaussian();
    void Init();

    void get_mean(double * _mean);
    void set_mean(double const* _mean);

    void get_covariance(double * _covariance);
    void set_covariance(double const* _covariance);

    void update_mean(double const* _mean);
    void update_covariance(double const* _covariance);

    double pdf(double const* x);
    double log_pdf(double const* x);
    REAL marginal_cdf(double const xcoord, size_t const marg_dim) const;
/*     REAL marginal_truncated_cdf(double const xcoord, size_t const marg_dim); */
    REAL inv_marginal_cdf(double const p, size_t const marg_dim) const;
/*     REAL inv_marginal_truncated_cdf(double const p, size_t const marg_dim); */
    void sample(double * x) const;
    void sample_conditioned(double const* x, double * y) { }
    void print_params(FILE * fh) const;
};


class ErrorEstimate;

class Posterior : public Integrand
{

public:

    ErrorEstimate * ee;
    double log_scaling_factor;

    double mode_point[4];
    bool zero_boundary[4];

    Posterior(ErrorEstimate * _ee,
              bool _may_underflow, size_t _ndim) :
        Integrand(_ndim, _may_underflow), ee(_ee)
    {
    }
    
    void initialize(double mode_tolerance, size_t max_function_evals, 
                    double const* initial_point, bool verbose);

    REAL pdf(double const* x);
    REAL log_pdf(double const* x);
    ErrorEstimate * model()
    {
        return ee;
    }

};


#endif // _INTEGRANDS_H
