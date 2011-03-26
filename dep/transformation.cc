#include "transformation.h"
#include "stats_tools.h"
#include "error_estimate.h"


#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

namespace Transformation
{


    double const sigmoid_truncation_epsilon = 1e-10;
    double const sigmoid_truncation_epsilon_inv = 
        1.0 - sigmoid_truncation_epsilon;

    /*
      filipe's iterative transformation: take a pair of additive inverses in one variable, and
      iteratively add another variable by multiplying by the sum of a pair:
      P(x) + Q(x) = 1
      P(x) + (P(y)+Q(y))Q(x) = 1
      (P(z) + Q(z))P(x) + (P(y)+Q(y))Q(x) = 1

      P(z)P(x) + Q(z)P(x) + P(y)Q(x) + Q(y)Q(x) = 1

      Let 
      P(x) + 1/(1 + exp(-x)), 
      Q(x) = exp(-x) / (1 + exp(-x))

      
    */
    void sigmoid_value_and_gradient(double const x[3],
                                    Transformation::SigmoidVals sg[3])
    {
        double enegx[3];

        double const epsilon = sigmoid_truncation_epsilon;
        double const epsilon_inv = 1.0 - epsilon;

        //here we assert that if the value p[d] is near 0 or 1,
        //the partial derivatives p'[d] and q'[d] we will deem
        //them to be zero, as if we'd already reached infinity
        //analytically, the sg and dirichlet cancel, and the
        //partial approaches a positive constant towards infinity.
        //but, for purposes of handling the dirichlet prior, we
        //make this correction.
        
        for (size_t d = 0; d != 3; ++d)
        {
            gsl_error_handler_t * original_handler = gsl_set_error_handler_off();
            enegx[d] = gsl_sf_exp(- x[d]);

            gsl_set_error_handler(original_handler);

            sg[d].p = 1.0 / (1.0 + enegx[d]);

            if (sg[d].p < epsilon)
            {
                sg[d].p = epsilon;
                sg[d].pgrad = 0.0;
                sg[d].pgrad_over_p = 0.0;
                sg[d].log_p = gsl_sf_log(epsilon);

                sg[d].q = epsilon_inv;
                sg[d].qgrad = 0.0;
                sg[d].qgrad_over_q = 0.0;
                sg[d].log_q = gsl_sf_log(epsilon_inv);
            }
            else if (sg[d].p > epsilon_inv)
            {
                sg[d].p = epsilon_inv;
                sg[d].pgrad = 0.0;
                sg[d].pgrad_over_p = 0.0;
                sg[d].log_p = gsl_sf_log(epsilon_inv);

                sg[d].q = epsilon;
                sg[d].qgrad = 0.0;
                sg[d].qgrad_over_q = 0.0;
                sg[d].log_q = gsl_sf_log(epsilon);
            }
            else 
            {
                sg[d].q = 1.0 - sg[d].p;
                sg[d].pgrad = (gsl_isinf(enegx[d]) == 1) ? 0.0 : enegx[d] / gsl_pow_2(1.0 + enegx[d]);
                sg[d].qgrad = - sg[d].pgrad;
                
                double log_1plus_enegx = gsl_sf_log_1plusx(enegx[d]);
                sg[d].log_p = -log_1plus_enegx;
                sg[d].log_q = - x[d] - log_1plus_enegx;
                
                sg[d].pgrad_over_p = sg[d].q;
                sg[d].qgrad_over_q = - sg[d].p;
            }

            assert(! isnan(sg[d].p));
            assert(! isnan(sg[d].q));
            assert(! isnan(sg[d].pgrad_over_p));
            assert(! isnan(sg[d].qgrad_over_q));
            assert(! isnan(sg[d].pgrad));
            assert(! isnan(sg[d].qgrad));
            assert(! isnan(sg[d].log_p));
            assert(! isnan(sg[d].log_q));
            
        }

    }                     

    void sigmoid_composition(Transformation::SigmoidVals const sg[3],
                             double * comp)
    {
        comp[0] = sg[2].p * sg[0].p;
        comp[1] = sg[2].q * sg[0].p;
        comp[2] = sg[1].p * sg[0].q;
        comp[3] = sg[1].q * sg[0].q;

        if (! (normalized(comp, 4, 1e-10) && all_positive(comp, 4)))
        {
            fprintf(stderr, "r3_to_composition_sigmoid: invalid composition\n");
            exit(2);
        }
        
    }

    //if the transformation has reached a point where the gradients are flat,
    //this is deemed to be 'infinity', which means the point represents the boundary
    void boundary_point(Transformation::SigmoidVals const sg[3],
                        bool * on_zero_boundary)
    {
        double const& e = Transformation::sigmoid_truncation_epsilon;

        on_zero_boundary[0] = sg[2].p == e || sg[0].p == e;
        on_zero_boundary[1] = sg[2].q == e || sg[0].p == e;
        on_zero_boundary[2] = sg[1].p == e || sg[0].q == e;
        on_zero_boundary[3] = sg[1].q == e || sg[0].q == e;

    }

    void sigmoid_gradient(Transformation::SigmoidVals const sg[3],
                          double comp_gradient[4][3])
    {
        comp_gradient[0][0] = sg[2].p * sg[0].pgrad;
        comp_gradient[0][1] = 0.0;
        comp_gradient[0][2] = sg[0].p * sg[2].pgrad;
        
        comp_gradient[1][0] = sg[2].q * sg[0].pgrad;
        comp_gradient[1][1] = 0.0;
        comp_gradient[1][2] = sg[0].p * sg[2].qgrad;

        comp_gradient[2][0] = sg[1].p * sg[0].qgrad;
        comp_gradient[2][1] = sg[0].q * sg[1].pgrad;
        comp_gradient[2][2] = 0.0;

        comp_gradient[3][0] = sg[1].q * sg[0].qgrad;
        comp_gradient[3][1] = sg[0].q * sg[1].qgrad;
        comp_gradient[3][2] = 0.0;

    }


    //compute log(~Dir(sigmoid(x)))
    double sigmoid_log_dirichlet(double const alpha[4],
                                 Transformation::SigmoidVals const sg[3])
    {
        double i[4];
        for (size_t c = 0; c != 4; ++c)
        {
            i[c] = alpha[c] - 1.0;
        }
        double value =
            + i[0] * (sg[2].log_p + sg[0].log_p)
            + i[1] * (sg[2].log_q + sg[0].log_p)
            + i[2] * (sg[1].log_p + sg[0].log_q)
            + i[3] * (sg[1].log_q + sg[0].log_q);

        return value;
    }
    
    //compute d/dx of log(~Dir(sigmoid(x)))
    void sigmoid_log_dirichlet_gradient(double const alpha[4],
                                        Transformation::SigmoidVals const sg[3],
                                        double * gradient)
    {
        double i[4];
        for (size_t c = 0; c != 4; ++c)
        {
            i[c] = alpha[c] - 1.0;
        }
        gradient[0] = 
            + i[0] * sg[0].pgrad_over_p 
            + i[1] * sg[0].pgrad_over_p 
            + i[2] * sg[0].qgrad_over_q 
            + i[3] * sg[0].qgrad_over_q;
        
        gradient[1] = i[2] * sg[1].pgrad_over_p + i[3] * sg[1].qgrad_over_q; 
        gradient[2] = i[0] * sg[2].pgrad_over_p + i[1] * sg[2].qgrad_over_q;
    }

    /*
      Let 

      Pinv(x) = -log(1/y - 1)

      c0 = P(z)P(x)
      c1 = Q(z)P(x)
      c2 = P(y)Q(x)
      c3 = Q(y)Q(x)

      so:

      x = Pinv(c0 + c1)
      y = Pinv(c2 / Q(x)) = Pinv(c2 / (c2 + c3)) = -log(1 + c3/c2 - 1)
      z = Pinv(c0 / P(x)) = Pinv(c0 / (c0 + c1)) = -log(1 + c1/c0 - 1)
      
    */

    //transform normalized 4D coordinates to r3 using a sigmoid
    //this needs to apply the truncation logic 
    //gsl_sf_log returns -nan if zero or close to zero with status 1,
    //inf if argument is larger than representable, with status 0
    void composition_to_r3_sigmoid(double const* c, double * r)
    {
        
        gsl_error_handler_t * original_handler = gsl_set_error_handler_off();
        double arg[3];
        arg[0] = 1.0 / (c[0] + c[1]) - 1.0;
        arg[1] = (c[3] == 0 && c[2] == 0) ? 1 : c[3] / c[2];
        arg[2] = (c[1] == 0 && c[0] == 0) ? 1 : c[1] / c[0];

        int status;
        gsl_sf_result result;
        status = gsl_sf_log_e(arg[0], & result);
        r[0] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        status = gsl_sf_log_e(arg[1], & result);
        r[1] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        status = gsl_sf_log_e(arg[2], & result);
        r[2] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        gsl_set_error_handler(original_handler);
    }


    double log_dirichlet(double const alpha[4],
                         double const x[4])
    {
        return
            + (alpha[0] - 1.0) * gsl_sf_log(x[0])
            + (alpha[1] - 1.0) * gsl_sf_log(x[1])
            + (alpha[2] - 1.0) * gsl_sf_log(x[2])
            + (alpha[3] - 1.0) * gsl_sf_log(x[3]);
    }


    

    void log_neg_posterior_aux(gsl_vector const* r, 
                               int eval_type, 
                               double * neg_log_value,
                               gsl_vector * neg_gradient_vec,
                               void * params)
    {
        double comp[4];
        double comp_gradient[4][3];
        double prior_gradient[3];
        double likelihood_gradient[4];
        double gradient[3];

        double log_prior;
        double log_likelihood;

        double x[3];
        x[0] = gsl_vector_get(r, 0);
        x[1] = gsl_vector_get(r, 1);
        x[2] = gsl_vector_get(r, 2);

        Transformation::PassingParams * pp = 
            static_cast<Transformation::PassingParams *>(params);

        Transformation::SigmoidVals sigmoid_vals[3];

        //smallest allowed distance from 0 or 1.  any closer to this,
        //and the sigmoid is considered to be equal to 0 or 1 with
        //derivative equal to 0

        sigmoid_value_and_gradient(x, sigmoid_vals);

        sigmoid_composition(sigmoid_vals, comp);

        //if (eval_type & Transformation::VALUE)
        if (1)
        {
            log_likelihood = pp->error_estimate->log_likelihood(comp) - pp->current_mode;
            log_prior = 
                sigmoid_log_dirichlet(pp->error_estimate->composition_prior_alphas,
                                      sigmoid_vals);

            *neg_log_value = -1.0 * (log_likelihood + log_prior);
        }

        //if (eval_type & Transformation::GRADIENT)
        if (1)
        {
            sigmoid_gradient(sigmoid_vals, comp_gradient);
            sigmoid_log_dirichlet_gradient(pp->error_estimate->composition_prior_alphas,
                                           sigmoid_vals, prior_gradient);
            
            pp->error_estimate->log_likelihood_gradient(comp, likelihood_gradient);

            for (size_t ri = 0; ri != 3; ++ri)
            {
                gradient[ri] = 0.0;
                for (size_t c = 0; c != 4; ++c)
                {
                    gradient[ri] += likelihood_gradient[c] * comp_gradient[c][ri];
                }
                gradient[ri] += prior_gradient[ri];
                gsl_vector_set(neg_gradient_vec, ri, gradient[ri] * -1.0);
            }
        }
        //if (1)
        if (0)
        {
            printf("etype: %i, "
                   "val: %10.8f, ll: %10.8f, prior: %10.8f, "
                   "r3: %10.8f, %10.8f, %10.8f, "
                   "prg: %10.8f, %10.8f, %10.8f, "
                   "grd: %10.8f, %10.8f, %10.8f, "
                   "composition: %10.8f, %10.8f, %10.8f, %10.8f\n",
                   eval_type,
                   *neg_log_value,
                   log_likelihood, log_prior,
                   x[0], x[1], x[2],
                   prior_gradient[0], prior_gradient[1], prior_gradient[2],
                   gradient[0], gradient[1], gradient[2],
                   comp[0], comp[1], comp[2], comp[3]);
        }
    }

    double log_neg_posterior_value(gsl_vector const* r, 
                                   void * params)
    {
        double neg_log_value;
        gsl_vector *neg_log_gradient = gsl_vector_alloc(3);

        log_neg_posterior_aux(r, Transformation::VALUE, &neg_log_value, neg_log_gradient, params);

        gsl_vector_free(neg_log_gradient);

        return neg_log_value;
    }


    void log_neg_posterior_gradient(const gsl_vector * r, 
                                    void * params, 
                                    gsl_vector * neg_gradient)
    {
        double neg_log_value;
        log_neg_posterior_aux(r, Transformation::GRADIENT, &neg_log_value, neg_gradient, params);
    }

    void log_neg_posterior_value_and_gradient(const gsl_vector * r, 
                                              void * params, 
                                              double * neg_log_value, 
                                              gsl_vector * neg_gradient)
    {
        log_neg_posterior_aux(r, Transformation::VALUE | Transformation::GRADIENT, 
                              neg_log_value, neg_gradient, params);
    }


    //     void log_neg_gradient_numerical(const gsl_vector * r, 
    //                                      void * params,
    //                                      double epsilon,
    //                                      gsl_vector * numerical_gradient)
    //     {
    //         gsl_vector * rdel = gsl_vector_alloc(3);
    //         double low, high;
        
    //         for (size_t d = 0; d != 3; ++d)
    //         {
    //             gsl_vector_memcpy(rdel, r);

    //             gsl_vector_set(rdel, d, gsl_vector_get(r, d) - epsilon);
    //             low = log_neg_posterior(rdel, params);

    //             gsl_vector_set(rdel, d, gsl_vector_get(r, d) + epsilon);
    //             high = log_neg_posterior(rdel, params);

    //             gsl_vector_set(numerical_gradient, d, (high - low) / (2.0 * epsilon));
            
    //         }
    //         gsl_vector_free(rdel);

    //     }

} // END namespace Transformation
