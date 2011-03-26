#include "integrands.h"

class Dirichlet : public SamplingFunction
{
    gsl_rng * seed;
    double * alpha;
    double Zlog2;
    double alpha0; // sum of alphas

public:
    Dirichlet(size_t _ndim, bool _may_underflow);
    ~Dirichlet();

    void update(double const* _alpha);
    void set_alpha0(double alpha0);
    double get_alpha0() { return this->alpha0; }

    void set_alphas_from_mode(double const* mode);
    void set_alphas_from_mean(double const* mean);

    void set_alphas_from_mode_or_bound(double const* mode,
                                       double const* alpha_lower_bound,
                                       bool const* is_zero_boundary);

    void set_alphas_from_mean_or_bound(double const* mean,
                                       double const* lower_bound);

    void lower_bound_alphas(double const* lower_bound);

    double pdf(double const* x);
    double log_pdf(double const* x);
    void sample(double * x) const;
    void sample_conditioned(double const* x_tau, 
                            double * x_star);
};


double
ran_dirichlet_lnpdf(const size_t K,
                    const double alpha[], const double theta[]);
