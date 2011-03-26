#include "metropolis.h"

#include <ctime>
#include <utility>
#include <deque>
#include <numeric>
#include <functional>
#include <algorithm>


#include <gmp.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>


#include "sampling.h"
#include "stats_tools.h"


//propose a new point z_star from a starting point z_tau
//z_tau is unused here since Metropolis-Hastings is using one
//unchanging proposal distribution
void Metropolis::propose(double const* z_tau, double * z_star)
{
    if (this->is_independence_chain_mh)
    {
        this->proposal->sample(z_star);
    }
    else
    {
        this->proposal->sample_conditioned(z_tau, z_star);
    }
}


double Metropolis::accept(double const* z_star, double const* z_tau)
{
    double ratio;

    if (this->integrand->may_underflow ||
        this->proposal->may_underflow)
    {
        double y_tau_log = this->integrand->log_pdf(z_tau);
        double y_star_log = this->integrand->log_pdf(z_star);
        double proposal_y_tau_log = this->proposal->log_pdf(z_tau);
        double proposal_y_star_log = this->proposal->log_pdf(z_star);

        assert(! isnan(y_tau_log));

        assert(! isnan(proposal_y_star_log));
        
        double numerator_log = y_star_log - proposal_y_star_log;
        double denominator_log = y_tau_log - proposal_y_tau_log;
        
        //master formula for metropolis hastings
        double ratio_log = numerator_log - denominator_log;
        
        assert(! isnan(ratio_log));
        
        //may be infinite.  this is okay.
        ratio = exp2f(ratio_log);
    }
    else
    {
        double y_tau = this->integrand->pdf(z_tau);
        double y_star = this->integrand->pdf(z_star);
        double proposal_y_tau = this->proposal->pdf(z_tau);
        double proposal_y_star = this->proposal->pdf(z_star);

        assert(y_tau != 0.0);
        assert(proposal_y_star != 0.0);

        double numerator = y_star / proposal_y_star;
        double denominator = y_tau / proposal_y_tau;
    
        //master formula for metropolis hastings
        ratio = numerator / denominator;

        assert(! isnan(ratio));
    }

    return ratio;

}



double * Metropolis::get_samples() const
{
    return sample_points;
}


void Metropolis::set_current_point(double const* point)
{
    std::copy(point, point + this->ndim, current_point);
}


double * Metropolis::get_current_point() const
{
    return this->current_point;
}


//iteratively refine the proposal gaussian by updating its mean and
//covariance according to sampled points.
/*
void 
Metropolis::tune_proposal(size_t num_points,
                          double max_tuning_iterations,
                          size_t averaging_window_size,
                          size_t front_back_window_distance,
                          size_t autocor_max_offset)
{

    //calculate mean, covariance, and autocorrelation
    size_t num_covariances = this->ndim * this->ndim;
    
    double * sample_mean = new double[this->ndim];
    double * estimated_mean = new double[this->ndim];
    double * sample_covar = new double[num_covariances];
    double * estimated_covar = new double[num_covariances];
    double * x = new double[this->ndim];
    double * x_next = new double[this->ndim];

    size_t autocor_buf_offset = front_back_window_distance;
    size_t autocor_bufsize = max_tuning_iterations + autocor_buf_offset;
    double * autocor_buffer = new double[autocor_bufsize];

    //2.0 is just any value above 1 that signals the non-convergence of the algorithm
    std::fill(autocor_buffer, autocor_buffer + autocor_buf_offset, 4.0);

    double * autocor_measure = autocor_buffer + autocor_buf_offset;

    double * starting_point = new double[this->ndim];
    std::copy(this->get_current_point(), this->get_current_point() + this->ndim,
              starting_point);

    if (this->integrand->pdf(starting_point) == 0)
    {
        fprintf(stderr, "Metropolis::tune_proposal: Error: starting point"
                " must yield a non-zero integrand value");
        exit(1);
    }

    this->proposal->update_mean(starting_point);
    this->proposal->get_mean(estimated_mean);
    this->proposal->get_covariance(estimated_covar);

    double proposal_avg;
    double front_mean, back_mean, front_stdev, back_stdev;

    for (size_t iteration = 0; iteration != max_tuning_iterations; ++iteration)
    {
        this->proposal->update_mean(estimated_mean);
        this->proposal->update_covariance(estimated_covar);
        
        //start the random walk from the mode point each time.
        std::copy(starting_point, starting_point + this->ndim, x);
        proposal_avg = 0.0;
        
        for (size_t n = 0; n != num_points; ++n)
        {
            double proposal_value = std::min(this->step(x, x_next), 1.0);

            proposal_avg = add_sample_to_mean(proposal_avg, n, proposal_value);

            std::copy(x_next, x_next + this->ndim, x);

            std::copy(x_next, x_next + this->ndim, 
                      this->get_samples() + (n * this->ndim));
        }
        
        fprintf(stdout, "iteration %Zu\n", iteration);

        autocor_measure[iteration] =
            total_autocorrelation_measure(this->get_samples(), this->ndim,
                                          num_points, autocor_max_offset);
        
        fprintf(stdout, "measure: %g, proposal: %g\n", 
                autocor_measure[iteration], proposal_avg);

        print_mean_covariance(stdout, estimated_mean, estimated_covar, this->ndim);

        fprintf(stdout, "\n\n");


        double * front_position = autocor_measure + iteration - averaging_window_size;
        double * back_position = front_position - averaging_window_size;

        front_mean = 
            gsl_stats_mean(front_position, 1, averaging_window_size);

        front_stdev = 
            gsl_stats_sd_m(front_position, 1, averaging_window_size, front_mean);

        back_mean = 
            gsl_stats_mean(back_position, 1, averaging_window_size);

        back_stdev = 
            gsl_stats_sd_m(back_position, 1, averaging_window_size, back_mean);

        if (fabs(front_mean - back_mean) < 0.01 * std::min(front_stdev, back_stdev)
            && front_mean < 1.0)
        {
            break;
        }
        
        
        multivariate_mean_covariance(this->get_samples(), this->ndim, 
                                     num_points, sample_mean, 
                                     sample_covar);


        double self_weight = 0.9;

        //shrink the covariance, leave the mean the same
        for (size_t d = 0; d != num_covariances; ++d)
        {
            estimated_covar[d] = 
                self_weight * estimated_covar[d] 
                + (1.0 - self_weight) * sample_covar[d];
        }
        for (size_t d = 0; d != this->ndim; ++d)
        {
            estimated_mean[d] = 
                self_weight * estimated_mean[d] 
                + (1.0 - self_weight) * sample_mean[d] ;
        }

    }

    delete sample_mean;
    delete estimated_mean;
    delete sample_covar;
    delete estimated_covar;
    delete x;
    delete x_next;
    delete starting_point;
    delete autocor_buffer;
}
*/





Metropolis::Metropolis(Integrand * _integrand,
                       SamplingFunction * _proposal,
                       size_t _ndim, 
                       bool _is_independence_chain,
                       size_t _num_points) 
    : integrand(_integrand), 
      proposal(_proposal), 
      ndim(_ndim), 
      is_independence_chain_mh(_is_independence_chain),
      num_points(_num_points)
{ 
    current_point = new double[ndim];
    sample_points = new double[ndim * _num_points];
}


Metropolis::~Metropolis()
{
    delete current_point;
    delete sample_points;
}


//do we need to record the y values here?
double Metropolis::step(double const* z_tau, double * z_tau_next)
{

    //update the markov chain at each step.
    double proposal_value;

    double * z_star = new double[this->ndim];
    double uniform;

    //from z_tau, propose a z_star
    this->propose(z_tau, z_star);

    assert(normalized(z_star, this->ndim, 1e-10));

    //from the integrand values, calculate and acceptance score
    proposal_value = this->accept(z_star, z_tau);

    uniform = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
            
    if (proposal_value > uniform)
    {
        //new jump is accepted; assign
        std::copy(z_star, z_star + this->ndim, z_tau_next);
    }

    else
    {
        //old position kept; assign
        std::copy(z_tau, z_tau + this->ndim, z_tau_next);
    }

    delete z_star;

    return proposal_value;

}


/* sample from the integrand using the Metropolis algorithm.
   populates the sample buffer */
void 
Metropolis::sample(size_t num_samples_to_take,
                   size_t burn_in,
                   size_t every_nth,
                   double * proposal_mean,
                   double * proposal_variance)
{
    if (num_samples_to_take > this->num_points)
    {
        fprintf(stderr, "Metropolis::sample:  %Zu samples requested exceeds "
                "%Zu samples allocated for this Metropolis instance"
                "Please call the constructor with larger number of samples", 
                num_samples_to_take, this->num_points);
        exit(1);
    }

//     if (integrand->pdf(this->current_point, this->ndim) == 0)
//     {
//         fprintf(stderr, "Metropolis::sample: Error: starting point"
//                 " must yield a non-zero integrand value");
//         exit(1);
//     }


    //update the markov chain at each step.
    srand(time(NULL));
    double proposal_value, proposal_value_truncated;
    size_t sample_count = 0;
    size_t step_count = 0;

    double * z_tau = new double[this->ndim];
    double * z_star = new double[this->ndim];

    std::copy(this->current_point, this->current_point + this->ndim, z_tau);
    std::pair<WEIGHTED_SAMPLE_MAP::iterator, bool> sample_insertion;
    std::vector<double> proposal_values;
    
    
    while (sample_count < num_samples_to_take)
    {

        proposal_value = this->step(z_tau, z_star);
        proposal_value_truncated = isinf(proposal_value) ? DBL_MAX : proposal_value;

        proposal_values.push_back(std::min(1.0, proposal_value_truncated));

        if (step_count > burn_in && step_count % every_nth == 0)
        {
            std::copy(z_tau, z_tau + this->ndim, this->get_samples() 
                      + (sample_count * this->ndim));
            ++sample_count;
        }

        ++step_count;

        //update current position
        std::copy(z_star, z_star + this->ndim, z_tau);

    }

    *proposal_mean = 
        gsl_stats_mean(&proposal_values[0], 1, proposal_values.size());

    *proposal_variance = 
        gsl_stats_variance_m(&proposal_values[0], 1, 
                             proposal_values.size(), *proposal_mean);

    this->set_current_point(z_tau);
 
    delete z_tau;
    delete z_star;
}


/*
other support functions
*/

void print_mean_covariance(FILE * fh, double const* mean, 
                           double const* covariance, 
                           size_t ndim)
{
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        fprintf(fh, "%20.18g\t", mean[d1]);
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            fprintf(fh, "\t%20.18g", covariance[d1 * ndim + d2]);
        }
        fprintf(fh, "\n");
    }
}





/*
void
Metropolis::fine_tune_proposal(Integrand * integrand,
                               ProposalDistribution * proposal,
                               size_t window_size,
                               size_t max_markov_steps,
                               double autocor_min_decrease,
                               size_t autocor_max_offset,
                               size_t initial_best_offset)
{
    if (window_size > this->num_points)
    {
        fprintf(stderr, "Metropolis::fine_tune_proposal: error: window size exceeds"
                " total number of samples for this Metropolis instance.");
        exit(1);
    }


    //fine tuning.  
    size_t num_sample_coords = this->ndim * window_size;

    double * sample_points = this->get_samples();
    double * mean = new double[this->ndim];
    double * covariance = new double[this->ndim * this->ndim];

    double * x = new double[this->ndim];
    double * x_next = new double[this->ndim];
    std::copy(this->current_point, this->current_point + this->ndim, x);

//     size_t current_best_offset = initial_best_offset;

    //automatically collect every nth sample, where n is determined
    //each time to be the best autocorrelation offset
    for (size_t n = 0; n != window_size; ++n)
    {
        this->step(integrand, proposal, x, x_next);
        std::copy(x_next, x_next + this->ndim, x);
        std::copy(x_next, x_next + this->ndim, 
                  sample_points + (n * this->ndim));
    }

    //calculate mean and covariance
    multivariate_mean_covariance(sample_points, this->ndim,
                                 window_size, mean, covariance);
    
    //initialize the sample window
    //std::deque<std::vector<double> > sample_window;
    std::vector<double> sample_window(num_sample_coords);
    std::copy(sample_points,
              sample_points + num_sample_coords,
              sample_window.begin());

    //this history 
    size_t autocor_history_offset = window_size;
    std::deque<double> autocor_history;
    autocor_history.resize(autocor_history_offset, 1.0);

    double autocor_decrease, autocor_measure;
    size_t step_number = 0;
    size_t current_window_size;

    do
    {
        //1.  perform one Metropolis-Hastings step
        
        this->step(integrand, proposal, x, x_next);
        
        //2.  update mean and covarance, and autocorrelation
        current_window_size = sample_window.size() / this->ndim;

        add_to_mean_covariance_matrix(x_next, this->ndim, 
                                      current_window_size, mean, covariance);
        
        sample_window.insert(sample_window.end(), x_next, x_next + this->ndim);
        current_window_size = sample_window.size() / this->ndim;

        remove_from_mean_covariance_matrix(&sample_window[0], this->ndim,
                                           current_window_size, mean, covariance);

        sample_window.erase(sample_window.begin(), sample_window.begin() + this->ndim);
        current_window_size = sample_window.size() / this->ndim;


        //3.  update proposal with mean and covariance from step 2.
        proposal->update_mean(mean);
        proposal->update_covariance(covariance);
        
        if (step_number % 1 == 0)
        {
            fprintf(stdout, "iteration %Zu\n", step_number);
            fprintf(stdout, "added point:");
            for (size_t d = 0; d != this->ndim; ++d)
            {
                fprintf(stdout, "\t%g", x_next[d]);
            }
            fprintf(stdout, "\n\n");
            
            print_mean_covariance(stdout, mean, covariance, this->ndim);
        }
        ++step_number;

        //4.  update markov chain
        std::copy(x_next, x_next + this->ndim, x);

        //5.  update autocorrelation history
        autocor_measure =
            total_autocorrelation_measure(&sample_window[0], this->ndim,
                                          current_window_size,
                                          autocor_max_offset);

        autocor_history.pop_front();
        autocor_history.push_back(autocor_measure);
        
        //6.  compute continuance criterion
        autocor_decrease = autocor_history.front() - autocor_history.back();
    }
    while ((autocor_decrease > autocor_min_decrease || autocor_measure > 1.0) &&
           step_number < max_markov_steps);

    delete mean;
    delete covariance;
    delete x;
    delete x_next;

}
*/
