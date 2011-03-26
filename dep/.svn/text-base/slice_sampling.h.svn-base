#ifndef _SLICE_SAMPLING_H
#define _SLICE_SAMPLING_H

#include <stdint.h>

#include <set>
#include <utility>

#include <gmp.h>

#include "tools.h"
#include "sampling.h"

/*
  Slice sampling as formulated by Radford M Neal and implemented here:
Introduction: Slice sampling is a Markov Chain Monte Carlo sampling
method that generates serially dependent set of points whose limiting
distribution is a sample from the function in question.  It samples an
N-dimensional function by traversing the N+1 dimensional space under
the function uniformly.  To do this, it starts at a given x, samples
from the marginal 

In Neal's words:

We sample alternately from the conditional distribution for y given
the current x, which is uniform over the interval (0, f (x)), and from
the conditional distribution for x given the current y, which is
uniform over the region S = {x :y <f(x)}, which I call the -Y´slice¡
defined by y.  Generating an independent point drawn uniformly from S
may still be difficult, in which case we can substitute some update
for x that leaves the uniform distribution over S invariant.


 */
class Integrand;


class SliceSampling
{
    size_t const ndim;
    size_t const nbits_per_dim;

    bool assume_log_integrand;
    size_t range_delta;

    size_t total_bits;

    mpz_t U, N, B, xprime, zero;
    mpf_t uniform;

    gmp_randstate_t rand_state;

    uint64_t * xcoord_grid;
    double * xcoord;
    int * coord_chunks;
    REAL yprime;


    bool initialized;

 public:

    SliceSampling(size_t const _ndim, size_t const nbits_per_dim, 
                  bool const _assume_log_integrand,
                  size_t _range_delta);

    ~SliceSampling();
    void Initialize();

    void grid_step(int const current_range, 
                   uint64_t const* xg,
                   uint64_t * xgprime);
    
    int step_in(Integrand * integrand,
                uint64_t const* xg,
                REAL y,
                int initial_range,
                uint64_t * xgp);

    int step_out(Integrand * integrand, 
                 uint64_t const* xg,
                 REAL const y,
                 int const initial_range);
    
    void initialize_starting_point(double const* starting_x, 
                                   size_t const ndim);
    

    REAL choose_auxiliary_coord(Integrand * integrand,
                                double const* x, 
                                size_t const ndim);

    void sample(Integrand * integrand,
                double const* starting_x,
                int const initial_range,
                size_t every_nth,
                double * sample_points_flat,
                size_t const num_samples);
    
};



//typedef std::binary_function<double *, size_t, REAL> Integrand;

//linearly scale ndim [0,1] (double) cube coords into (uint64_t)
//[0,max_64bit_int] grid_coords
void expand_to_grid_coords(double const* cube_coord, size_t const ndim, 
                           uint64_t * grid_coord, size_t const num_chunks);

//linearly scale ndim (uint64_t) [0,max_64bit_int] grid_coords into
//[0,1] (double) cube coords
void contract_from_grid_coords(uint64_t const* grid_coord, size_t const ndim,
                               double * cube_coord, size_t const num_chunks);


#endif // _SLICE_SAMPLING_H
