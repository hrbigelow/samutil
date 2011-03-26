#include <climits>

#include <gsl/gsl_sf_log.h>

#include "slice_sampling.h"

#include "tools.h"
#include "hilbert.h"
#include "sampling.h"
#include "integrands.h"

SliceSampling::SliceSampling(size_t const _ndim, 
                             size_t const _nbits_per_dim,
                             bool const _assume_log_integrand,
                             size_t _range_delta) : 
    ndim(_ndim),
    nbits_per_dim(_nbits_per_dim),
    assume_log_integrand(_assume_log_integrand),
    range_delta(_range_delta),
    initialized(false)
{
    xcoord_grid = new uint64_t[ndim];
    xcoord = new double[ndim];
    coord_chunks = new int[nbits_per_dim];
}


SliceSampling::~SliceSampling()
{
    if (xcoord != NULL)
    {
        delete xcoord;
        xcoord = NULL;
    }
    if (xcoord_grid != NULL)
    {
        delete xcoord_grid;
        xcoord_grid = NULL;
    }
    if (coord_chunks != NULL)
    {
        delete coord_chunks;
        coord_chunks = NULL;
    }
}


void SliceSampling::Initialize()
{
    this->total_bits = this->ndim * this->nbits_per_dim;

    mpz_init(this->U);
    mpz_init(this->N);
    mpz_init(this->B);
    mpz_init(this->zero);
    mpz_init(this->xprime);
    mpf_init(this->uniform);

    mpz_set_ui(this->B, 0);
    mpz_set_ui(this->zero, 0);
    mpz_set_d(this->B, pow(2,this->total_bits));

    gmp_randinit_default(this->rand_state);

    //size_t const bits = sizeof(double) * 8;
    this->initialized = true;
}
    



//linearly scale ndim [0,1] (double) cube coords into (uint64_t)
//[0,max_64bit_int] grid_coords
void expand_to_grid_coords(double const* cube_coord, size_t const ndim, 
                           uint64_t * grid_coord, size_t const nbits_per_dim)
{
    //rescale the number from [0,1] to [0, max_integer]
    double factor = static_cast<double>((static_cast<uint64_t>(1)<<(nbits_per_dim+1)) - 1);

    for (size_t i = 0; i != ndim; ++i)
    {
        grid_coord[i] = static_cast<uint64_t>(cube_coord[i] * factor);
    }
}


//linearly scale ndim (uint64_t) [0,max_64bit_int] grid_coords into
//[0,1] (double) cube coords
void contract_from_grid_coords(uint64_t const* grid_coord, size_t const ndim,
                               double * cube_coord, size_t const nbits_per_dim)
{
    double factor = static_cast<double>((static_cast<uint64_t>(1)<<(nbits_per_dim+1)) - 1);

    for (size_t i = 0; i != ndim; ++i)
    {
        cube_coord[i] = static_cast<double>(grid_coord[i]) / factor;
    }
}



/* sample a bits interval of current_range in the neighborhood of x,
   storing the result in xprime
 */
/*
void SliceSampling::sample_bits_interval(int const current_range, 
                                         mpz_t const& x,
                                         uint64_t * xg,
                                         mpz_t * xnew,
                                         uint64_t * xgnew)
{

    mpz_urandomb(this->N, this->rand_state, current_range);
        
    mpz_sub(*xnew, x, this->U);
        
    if (mpz_cmp(*xnew, this->zero) < 0)
    {
        //xprime < 0: wrap it
        mpz_add(*xnew, *xnew, this->B);
        int cmp = mpz_cmp(*xnew, this->zero);
        assert(cmp > 0);
    }
    int_to_Hilbert(*xnew, this->nbits_per_dim, xgnew, this->ndim);

    mpz_xor(*xnew, *xnew, this->N);
    mpz_add(*xnew, *xnew, this->U);
    if (mpz_cmp(*xnew, this->B) > 0)
    {
        //xprime > this->B:  wrap it
        mpz_sub(*xnew, *xnew, this->B);
        int cmp = mpz_cmp(*xnew, this->B);
        assert(cmp < 0);
    }

}
*/


//randomly step in the integer grid with a step size up to <current_range> bits
void SliceSampling::grid_step(int const current_range, 
                              uint64_t const* xg,
                              uint64_t * xgprime)
{

    mpz_t x;
    mpz_init(x);

    uint64_t * xgu = new uint64_t[this->ndim];
    uint64_t * xgu_prime = new uint64_t[this->ndim];

    Hilbert_to_int(xg, this->ndim, this->nbits_per_dim, x);

    mpz_urandomb(this->N, this->rand_state, current_range);
    mpz_sub(x, x, this->U);
    if (mpz_cmp(x, this->zero) < 0)
    {
        //xprime < 0: wrap it
        mpz_add(x, x, this->B);
        int cmp = mpz_cmp(x, this->zero);
        assert(cmp > 0);
    }
    int_to_Hilbert(x, this->nbits_per_dim, xgu, this->ndim);

    mpz_xor(x, x, this->N);
    int_to_Hilbert(x, this->nbits_per_dim, xgu_prime, this->ndim);

//     mpz_add(*xnew, *xnew, this->U);
//     if (mpz_cmp(*xnew, this->B) > 0)
//     {
//         //xprime > this->B:  wrap it
//         mpz_sub(*xnew, *xnew, this->B);
//         int cmp = mpz_cmp(*xnew, this->B);
//         assert(cmp < 0);
//     }

    for (size_t d = 0; d != this->ndim; ++d)
    {
        xgprime[d] = xg[d] + xgu_prime[d] - xgu[d];
    }

    delete xgu;
    delete xgu_prime;
}


/* sample from the search interval until a point is within the slice
   (integrand value is above y).  at each point not within the slice,
   shrink the interval so as to retain the initial point x within it.
   interval starts as a neighborhood of x with initial_range bits of
   width.  updates xcoord, xcoord_grid and yprime with latest values
*/
int SliceSampling::step_in(Integrand * integrand,
                           uint64_t const* xg,
                           REAL y,
                           int initial_range,
                           uint64_t * xgp)
{

    int current_range = initial_range + this->range_delta;
    mpz_urandomb(this->U, this->rand_state, this->total_bits);
    double * xrp = new double[this->ndim];
    double yp;

    do
    {
        current_range -= this->range_delta;
        this->grid_step(current_range, xg, xgp);
        contract_from_grid_coords(xgp, this->ndim, xrp, this->nbits_per_dim);
        if (this->assume_log_integrand)
        {
            yp = integrand->log_pdf(xrp);
        }
        else
        {
            yp = integrand->pdf(xrp);
        }
    }
    while (yp < y && (! std::equal(xg, xg + this->ndim, xgp)));
    
    delete xrp;

    return current_range;
}



/* expand the search interval until a test point is found that lies
   outside (below) the slice defined by level y.  start in the
   neighborhood of x with a neighborhood of intitial_range
   bits. returns the range found.
*/
int SliceSampling::step_out(Integrand * integrand, 
                            uint64_t const* xg,
                            REAL const y,
                            int const initial_range)
{

    int current_range = initial_range - this->range_delta;
    mpz_urandomb(this->U, this->rand_state, this->total_bits);
    uint64_t * xgp = new uint64_t[this->ndim];
    double * xrp = new double[this->ndim];
    double yp;

    do
    {
        current_range += this->range_delta;

        this->grid_step(current_range, xg, xgp);
        contract_from_grid_coords(xgp, this->ndim, xrp, this->nbits_per_dim);
        if (this->assume_log_integrand)
        {
            yp = integrand->log_pdf(xrp);
        }
        else
        {
            yp = integrand->pdf(xrp);
        }
    }
    while (current_range != static_cast<int>(this->total_bits) && yp >= y);

    delete xgp;
    delete xrp;

    return current_range;
}


/* initialize xcoord, xcoord_grid (the projection of xcoord onto an
   integer grid), and xprime (xcoord_grid mapped onto the hilbert
   curve) from starting_x (plain unit_hypercube x coordinates)
   
 */
void SliceSampling::initialize_starting_point(double const* starting_x, 
                                              size_t const ndim)
{

    mpf_urandomb(this->uniform, this->rand_state, 61);

    std::copy(starting_x, starting_x + ndim, this->xcoord);
    
    expand_to_grid_coords(this->xcoord, ndim, this->xcoord_grid, 
                          this->nbits_per_dim);

    Hilbert_to_int(this->xcoord_grid, ndim, this->nbits_per_dim, this->xprime);
    
}    


/* selects an auxiliary coordinate between U(0, f(x)).  If
   assume_log_integrand is true, the auxiliary coordinate is also
   sampled in log space, and the integrand given is assumed to be
   log(integrand-of-interest)
 */
REAL SliceSampling::choose_auxiliary_coord(Integrand * integrand,
                                           double const* x, 
                                           size_t const ndim)
{

    mpf_urandomb(this->uniform, this->rand_state, 64);
    REAL y;
    
    if (this->assume_log_integrand)
    {
        y = integrand->log_pdf(x) + gsl_sf_log(mpf_get_d(this->uniform));
    }
    else
    {
        y = integrand->pdf(x) * mpf_get_d(this->uniform);
    }
    return y;
}


/* xprime and yprime are the proposed new point for slice sampling they
   are only accepted if yprime < y.  once accepted, the slice sampling
   continues, using (xprime, yprime) as the new (x, y)
*/
void SliceSampling::sample(Integrand * integrand,
                           double const* starting_x,
                           int initial_range,
                           size_t every_nth,
                           double * sample_points_flat,
                           size_t num_samples)
{

    if (! this->initialized)
    {
        fprintf(stderr, "Must first initialize SliceSampling\n");
        exit(1);
    }

    if (static_cast<size_t>(initial_range) > this->total_bits)
    {
        fprintf(stderr, "initial_range (%i) must be <= total bits (%Zu)\n",
                initial_range, this->total_bits);
        exit(1);
    }

    uint64_t * xg = new uint64_t[this->ndim];
    uint64_t * xgp = new uint64_t[this->ndim];
    double * xrp = new double[this->ndim];

    //initialize
    int current_range = initial_range;
    std::copy(starting_x, starting_x + this->ndim, xrp);
    expand_to_grid_coords(xrp, ndim, xg, this->nbits_per_dim);
    REAL y = this->choose_auxiliary_coord(integrand, xrp, this->ndim);

    size_t sample_count = 0;

    for (size_t si = 0; si != num_samples * every_nth; ++si)
    {
        current_range = initial_range;
        
        current_range = this->step_out(integrand, xg, y, current_range);
//         fprintf(stdout, "out to %i", current_range);
        
        current_range = this->step_in(integrand, xg, y, current_range, xgp);
//         fprintf(stdout, ", in to %i\n", current_range);

        contract_from_grid_coords(xgp, this->ndim, xrp, this->nbits_per_dim);

        if (si % every_nth == 0)
        {
            std::copy(xrp, xrp + this->ndim, 
                      sample_points_flat + (sample_count * this->ndim));
            ++sample_count;
        }

        //update markov chain
        std::copy(xgp, xgp + this->ndim, xg);

        //choose y coordinate from new x coordinate
        y = this->choose_auxiliary_coord(integrand, xrp, this->ndim);
        
    }

    delete xg;
    delete xgp;
    delete xrp;

}
