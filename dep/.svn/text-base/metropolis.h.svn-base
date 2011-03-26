#ifndef _METROPOLIS_H
#define _METROPOLIS_H

#include <gsl/gsl_randist.h>

#include "sampling.h"
#include "integrands.h"

class Metropolis
{

    Integrand * integrand;
    SamplingFunction * proposal;

    size_t ndim;

    //if true, use 'independence chain MH sampling and thus produce
    //proposals independent of z_tau
    bool is_independence_chain_mh;
    double * current_point;
    double * sample_points;

 public:

    size_t num_points;

    Metropolis(Integrand * _integrand,
               SamplingFunction * _proposal,
               size_t _ndim, 
               bool _is_independence_chain,
               size_t total_sample_points);

    ~Metropolis();

    double * get_samples() const;
    void set_current_point(double const* point);
    double * get_current_point() const;

    void propose(double const* z_tau, double * z_star);
    double accept(double const* z_star, double const* z_tau);

    void sample(size_t const num_samples,
                size_t const burn_in,
                size_t const every_nth,
                double * proposal_mean,
                double * proposal_variance);
    
    double step(double const* z_tau, double * z_tau_next);
    
    void tune_proposal(size_t num_points,
                       double max_tuning_iterations,
                       size_t averaging_window_size,
                       size_t front_back_window_distance,
                       size_t autocor_max_offset);
    
};


void print_mean_covariance(FILE * fh, 
                           double const* mean, 
                           double const* covariance, 
                           size_t ndim);


#endif // _METROPOLIS_H
