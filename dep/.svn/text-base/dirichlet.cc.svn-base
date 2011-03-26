#include "dirichlet.h"
#include "stats_tools.h"
#include "transformation.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>

#include <functional>
#include <algorithm>
#include <numeric>

Dirichlet::Dirichlet(size_t _ndim, bool _may_underflow)
    : SamplingFunction(_ndim, _may_underflow)
{
    this->alpha = new double[this->ndim];
    this->seed = gsl_rng_alloc(gsl_rng_taus);
    
}

Dirichlet::~Dirichlet()
{
    delete this->alpha;
    gsl_rng_free(this->seed);
}


void Dirichlet::update(double const* _alpha)
{
    std::copy(_alpha, _alpha + this->ndim, this->alpha);
    this->alpha0 = 0;

    for (size_t d = 0; d != this->ndim; ++d)
    {
        this->alpha0 += this->alpha[d];

        if (this->alpha[d] <= 0.0)
        {
            fprintf(stderr, "Dirichlet::update: alphas should be all > 0\n");
            exit(1);
        }
    }
}


void Dirichlet::set_alpha0(double alpha0)
{
    this->alpha0 = alpha0;
}


void Dirichlet::set_alphas_from_mode(double const* mode)
{
    for (size_t d = 0; d != this->ndim; ++d)
    {
        this->alpha[d] = mode[d] * (this->alpha0 - this->ndim) + 1;
    }
}


/* 
   The purpose of this function is to adjust the proposal Dirichlet
   alphas to fit the posterior as best as possible.  Since the
   posterior combines the concave Dirichlet prior with highly variable
   evidence, the posterior may be convex in some dimensions (the
   dimensions with base observations) but still concave in others.

   Note that a posterior with no base observations reduces to a pure
   Dirichlet, which is the prior. In that case, simply knowing the 

   1) if the mode_or_peak is on the zero boundary, the alpha should be
   set to a lower-bound value.

   2) if the mode_or_peak is on the ones boundary, it means the
      ones-complement dimension

   2) otherwise it is concave in that dimension
 */
void Dirichlet::set_alphas_from_mode_or_bound(double const* mode_or_peak,
                                              double const* alpha_lower_bound,
                                              bool const* is_zero_boundary)
{
    size_t num_observed_dims = 0;

    assert(this->alpha0 > 1);

    for (size_t d = 0; d != this->ndim; ++d)
    {
        if (! is_zero_boundary[d])
        {
            num_observed_dims++;
        }
    }

    //assert(this->alpha0 > num_observed_dims);

    for (size_t d = 0; d != this->ndim; ++d)
    {
        if (is_zero_boundary[d])
        {
            this->alpha[d] = 0.0;
        }
        else
        {
            this->alpha[d] = 
                mode_or_peak[d] * (this->alpha0 - num_observed_dims) + 1;
        }
    }
    this->lower_bound_alphas(alpha_lower_bound);
}



void Dirichlet::set_alphas_from_mean(double const* mean)
{
    for (size_t d = 0; d != this->ndim; ++d)
    {
        this->alpha[d] = mean[d] * this->alpha0;
    }
}


void Dirichlet::set_alphas_from_mean_or_bound(double const* mean,
                                              double const* lower_bound)
{
    this->set_alphas_from_mean(mean);
    this->lower_bound_alphas(lower_bound);
}


//up-adjust alphas that are less than the lower_bound given
//re-adjust all other alphas, keeping ratios intact, so that
//the sum-of-alphas coincides with alpha0
void Dirichlet::lower_bound_alphas(double const* lower_bound)
{
    double qual_alpha0 = 0.0;
    double new_alpha0 = 0.0;

    //raise all unsufficient alphas
    for (size_t d = 0; d != this->ndim; ++d)
    {
        if (this->alpha[d] < lower_bound[d])
        {
            this->alpha[d] = lower_bound[d];
        }
        else
        {
            qual_alpha0 += this->alpha[d];
        }
        new_alpha0 += this->alpha[d];
    }
    double adjust = 1.0 + (this->alpha0 - new_alpha0) / qual_alpha0;
    for (size_t d = 0; d != this->ndim; ++d)
    {
        if (this->alpha[d] != lower_bound[d])
        {
            this->alpha[d] *= adjust;
        }
    }
    double sum = std::accumulate(this->alpha, this->alpha + 4, 0.0);
    int compare = gsl_fcmp(sum, this->alpha0, 1e-5);
    assert(compare == 0);
    
}


double Dirichlet::pdf(double const* x)
{
    return gsl_ran_dirichlet_pdf(this->ndim, this->alpha, x);
}



/*
  Reimplementation of GSL's function gsl_ran_dirichlet_lnpdf that handles
  cases where alpha_i == 1 or alpha_i == 0.
  When alpha_i == 1, the exponent is zero, allowing the opportunity for
  theta to be zero.  but, when alpha_i == 0, gsl_sf_lngamma is undefined.
  
 */
double
ran_dirichlet_lnpdf(const size_t K,
                    const double alpha[], const double theta[])
{
    /*We calculate the log of the pdf to minimize the possibility of overflow */
    size_t i;
    double log_p = 0.0;
    double sum_alpha = 0.0;

    for (i = 0; i < K; i++)
    {
        log_p += alpha[i] == 1.0 ? 0.0 : (alpha[i] - 1.0) * log(theta[i]);
    }
  
    for (i = 0; i < K; i++)
    {
        sum_alpha += alpha[i];
    }

    log_p += gsl_sf_lngamma (sum_alpha);

    for (i = 0; i < K; i++)
    {
        log_p -= gsl_sf_lngamma (alpha[i]);
    }

    return log_p;
}


double Dirichlet::log_pdf(double const* x)
{
    return Transformation::log_dirichlet(this->alpha, x);
}


void Dirichlet::sample(double * x_star) const
{
    gsl_ran_dirichlet(seed, this->ndim, this->alpha, x_star);
    double r3[3];
    Transformation::composition_to_r3_sigmoid(x_star, r3);

    Transformation::SigmoidVals sigmoid_vals[3];
    Transformation::sigmoid_value_and_gradient(r3, sigmoid_vals);
    
    Transformation::sigmoid_composition(sigmoid_vals, x_star);
}


void Dirichlet::sample_conditioned(double const* x_tau, 
                                   double * x_star)
{
    this->set_alphas_from_mode(x_tau);
    gsl_ran_dirichlet(seed, this->ndim, this->alpha, x_star);
}
