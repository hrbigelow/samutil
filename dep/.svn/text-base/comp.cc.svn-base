#include <vector>
#include <utility>
#include <numeric>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include <sys/timeb.h>

#include "tools.h"
#include "error_estimate.h"
#include "integrands.h"
#include "metropolis.h"
#include "sampling.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "dirichlet.h"
#include "slice_sampling.h"
#include "stats_tools.h"
#include "nucleotide_reader.h"
#include "nucleotide_stats.h"
#include "base_qual_strand_reader.h"
#include "simulation.h"
#include "anomaly_tools.h"
#include "transformation.h"


#include "usage_strings.h"

int comp_usage()
{
    fprintf(stderr, 
            "\nUsage: dep comp [options] jpd_params input.pileup output.comp [output.points]\n" 
            "Options:\n\n"
            "-d STRING   %s\n"
            "-l STRING   %s\n"
            "-t INT      number of sample points used for tuning proposal distribution [1000]\n"
            "-f INT      number of sample points used for final quantiles estimation [10000]\n"
            "-a INT      target autocorrelation offset.  once reached, proposal tuning halts [6]\n"
            "-i INT      maximum number of proposal tuning iterations [10]\n"
            "-s INT      number of loci simulations for estimating anomaly score (if zero, no estimate provided) [0]\n"
            "-m FLOAT    autocorrelation maximum value [6]\n"
            "-Q FILE     quantiles file [\"0.005 0.05 0.5 0.95 0.995\\n\"]\n"
            "-p FILE     alpha values for dirichlet prior [\"0.1 0.1 0.1 0.1\\n\"]\n"
            "-q INT      %s\n"
            "-v <empty>  if present, be verbose [absent]\n"
            "\n", Usage::data_adaptor, Usage::label_string, Usage::quality_string
            );
    return 1;
}

extern char *optarg;
extern int optind;



int main_comp(int argc, char ** argv)
{

    char input_data_type[100];
    strcpy(input_data_type, "base_qual_strand");

    char label_string[100];
    strcpy(label_string, "comp");

    size_t tuning_num_points = 1e3;
    size_t final_num_points = 1e4;

    size_t target_autocor_offset = 6;
    size_t max_tuning_iterations = 10;

    size_t nsim_loci = 0;

    double autocor_max_value = 6;
    char quantiles_file[100];
    strcpy(quantiles_file, "/dev/null");

    char prior_alphas_file[100];
    strcpy(prior_alphas_file, "/dev/null");


    size_t min_quality_score = 5;

    double default_prior_alpha = 0.1;

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* posterior_output_file;
    char const* cdfs_output_file;

    bool verbose = false;

    char c;
    while ((c = getopt(argc, argv, "d:l:t:f:a:s:i:m:Q:p:q:v")) >= 0)
    {
        switch(c)
        {
        case 'd': strcpy(input_data_type, optarg); break;
        case 'l': strcpy(label_string, optarg); break;
        case 't': tuning_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'f': final_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'a': target_autocor_offset = static_cast<size_t>(atof(optarg)); break;
        case 'i': max_tuning_iterations = static_cast<size_t>(atof(optarg)); break;
        case 's': nsim_loci = static_cast<size_t>(atof(optarg)); break;
        case 'm': autocor_max_value = atof(optarg); break;
        case 'Q': strcpy(quantiles_file, optarg); break;
        case 'p': strcpy(prior_alphas_file, optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'v': verbose = true; break;
        default: return comp_usage(); break;
        }
    }
    if (argc - optind != 3)
    {
        return comp_usage();
    }

    jpd_data_params_file = argv[optind];
    pileup_input_file = argv[optind + 1];
    posterior_output_file = argv[optind + 2];

    cdfs_output_file = (optind + 3 < argc) ? argv[optind + 3] : "/dev/null";
    
    //parse mass fractions file
    double * quantiles;
    size_t num_quantiles;

    if (strcmp(quantiles_file, "/dev/null") == 0)
    {
        quantiles = new double[5];
        quantiles[0] = 0.005;
        quantiles[1] = 0.05;
        quantiles[2] = 0.5;
        quantiles[3] = 0.95;
        quantiles[4] = 0.995;
        num_quantiles = 5;
    }
    else
    {
        quantiles = ParseNumbersFile(quantiles_file, & num_quantiles);
    }

    double * prior_alphas;
    size_t num_prior_alphas;

    if (strcmp(prior_alphas_file, "/dev/null") == 0)
    {
        prior_alphas = new double[4];
        std::fill(prior_alphas, prior_alphas + 4, default_prior_alpha);
    }
    else
    {
        prior_alphas = ParseNumbersFile(prior_alphas_file, & num_prior_alphas);
    }

    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);
    
    size_t full_ndim = 4;
    size_t truncated_ndim = 3;

    FILE * posterior_output_fh = fopen(posterior_output_file, "w");
    if (posterior_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open posterior_output_file %s\n",
                posterior_output_file);
    }

    FILE * cdfs_output_fh;
    if (strcmp(cdfs_output_file, "/dev/null") == 0)
    {
        cdfs_output_fh = NULL;
    }
    else
    {
        cdfs_output_fh = fopen(cdfs_output_file, "w");
        if (cdfs_output_fh == NULL)
        {
            fprintf(stderr, "Couldn't open cdfs_output_file %s\n",
                    cdfs_output_file);
        }
    }

    NucleotideReader * data_reader;

    if (strcmp(input_data_type, "base_qual_strand") == 0)
    {
        data_reader = new BaseQualStrandReader;
    }
    //insert other data handlers here
    else
    {
        fprintf(stderr,
                "Currently only 'base_qual_strand' data is supported\n");
        exit(5);
    }

    data_reader->initialize(pileup_input_file);

    //we are integrating the actual posterior
    char const* dimension_labels[] = { "A", "C", "G", "T" };
    char line_label[1000];

    double mode_tolerance = 1e-60;
    size_t max_modefinding_iterations = 3000;
    double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };


    bool use_independence_chain_mh = true;
    bool may_underflow = true;

    size_t effective_depth;
    
    double prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);

    ErrorEstimate model;
    model.set_composition_prior_alphas(prior_alphas);
    
    //1. Compile JPD_DATA from counts file
    JPD_DATA global_counts = parse_jpd_rdb_file(jpd_data_params_file);
    
    //2. initialize NucleotideStats from JPD_DATA
    NucleotideStats params(global_counts.size());
    params.initialize(global_counts);

    //3. Set model parameters
    model.set_model_params(&params);

    //4. Construct posterior
    Posterior posterior(&model, may_underflow, full_ndim);


    while (data_reader->more_loci())
    {
        //1. get entire locus data
        LocusSummary locus = 
            data_reader->get_next_locus(params, static_cast<void const*>(& min_quality_score));

        //divide locus data to plus and minus-strand data
        double pos_anomaly_score =
            strand_locus_anomaly_score(posterior, global_counts,
                                       locus, data_reader, '+', verbose);
        
        double neg_anomaly_score =
            strand_locus_anomaly_score(posterior, global_counts,
                                       locus, data_reader, '-', verbose);
        
        
        effective_depth = locus.read_depth;

        Dirichlet dirichlet(full_ndim, may_underflow);

        posterior.model()->set_locus_data(&locus);
        posterior.model()->get_model_params()->initialize(global_counts);

        posterior.initialize(mode_tolerance, max_modefinding_iterations, 
                             initial_point, verbose);



        Metropolis metropolis(&posterior, &dirichlet, full_ndim, 
                              use_independence_chain_mh, final_num_points);


        double autocor_max_offset = 6;

        metropolis.set_current_point(posterior.mode_point);
        
        double estimated_mean[full_ndim];

        size_t initial_autocor_offset = 30;

        size_t best_autocor_offset = initial_autocor_offset;

        dirichlet.set_alpha0(locus.read_depth + prior_alpha0);

        if (dirichlet.get_alpha0() > 1)
        {
            dirichlet.set_alphas_from_mode_or_bound
                (posterior.mode_point,
                 posterior.ee->composition_prior_alphas,
                 posterior.zero_boundary);
        }
        else
        {
            dirichlet.update(posterior.ee->composition_prior_alphas);
        }

        //metropolis hastings
        double proposal_mean, proposal_variance;

        size_t cumul_autocor_offset;
        for (size_t iter = 0; iter != max_tuning_iterations; ++iter)
        {
            cumul_autocor_offset = 1;
            size_t current_autocor_offset = 1;
            for (size_t i = 0; i != 3; ++i)
            {
                //sample more and more thinly, starting from every 1'th
                metropolis.sample(tuning_num_points, 0, cumul_autocor_offset,
                                  &proposal_mean, &proposal_variance);
                
                current_autocor_offset =
                    best_autocorrelation_offset(metropolis.get_samples(),
                                                full_ndim, tuning_num_points,
                                                autocor_max_offset, autocor_max_value);

                if (current_autocor_offset == 1)
                {
                    break;
                }

                cumul_autocor_offset *= current_autocor_offset;
                if (verbose)
                {
                    fprintf(stdout, "MH: current: %Zu, cumulative: %Zu, position: %Zu\n", 
                            current_autocor_offset, cumul_autocor_offset, 
                            locus.position);
                    fflush(stdout);
                }
            }

            //here it doesn't make sense to take all of the samples if
            //the best_autoocr offset isn't 1
            multivariate_mean(metropolis.get_samples(), full_ndim,
                              tuning_num_points, estimated_mean);

            if (verbose)
            {
                fprintf(stdout, "MH: %Zu, position: %Zu, proposal mean: %g, "
                        "proposal variance: %g", cumul_autocor_offset,
                        locus.position, proposal_mean, proposal_variance);
                fprintf(stdout, ", estimated mean:");
                for (size_t d = 0; d != full_ndim; ++d)
                {
                    fprintf(stdout, "\t%g", estimated_mean[d]);
                }
                fprintf(stdout, "\n");
                fflush(stdout);
            }

            if (cumul_autocor_offset <= target_autocor_offset)
            {
                break;
            }

            dirichlet.set_alphas_from_mean_or_bound(estimated_mean,
                                                    posterior.ee->composition_prior_alphas);

            
        }


        if (cumul_autocor_offset <= target_autocor_offset)
        {
            //metropolis hastings succeeded
            metropolis.sample(final_num_points, 0, cumul_autocor_offset,
                              &proposal_mean, &proposal_variance);

            sprintf(line_label, "%s\t%s\t%s\t%Zu\t%c\t%Zu\t%Zu\t%5.5lf\t%5.5lf", label_string, 
                    "MH", locus.reference, 
                    locus.position, 
                    locus.reference_base, 
                    locus.read_depth, effective_depth,
                    pos_anomaly_score,
                    neg_anomaly_score);

            double * sample_points_buf = new double[final_num_points * 4];
            std::vector<double *> sample_points_sortable(final_num_points);

            for (size_t i = 0; i != final_num_points; ++i)
            {
                sample_points_sortable[i] = metropolis.get_samples() + (i * full_ndim);
            }
            
//             //approximate marginal modes (expensive!)
//             std::fill(marginal_modes, marginal_modes + 4, -1.0);
//             if (marginal_mode_num_points > 0)
//             {
//                 for (size_t d = 0; d != full_ndim; ++d)
//                 {
//                     marginal_modes[d] = 
//                         window_averaged_mode(& sample_points_sortable, d, 
//                                              marginal_mode_num_points);
//                 }
//             }

            print_quantiles(posterior_output_fh, & sample_points_sortable, 
                            posterior.mode_point,
                            line_label, dimension_labels, "+", quantiles,
                            num_quantiles, full_ndim);

            if (cdfs_output_fh != NULL)
            {
                print_numerical_cdfs(cdfs_output_fh, line_label, & sample_points_sortable, full_ndim);
            }

            fflush(posterior_output_fh);
            delete sample_points_buf;

            continue;

        }

        if (verbose)
        {
            fprintf(stderr, 
                    "%s %Zu (depth %Zu) locus skipped"
                    ", Metropolis Hastings failed.\n",
                    locus.reference, locus.position, effective_depth);
        }

        continue;


        //metropolis hastings failed.  try slice sampling
        size_t num_bits_per_dim = 62;
        bool is_log_integrand = true;

        size_t initial_sampling_range = 62 * truncated_ndim;

        SliceSampling slice_sampler(truncated_ndim, num_bits_per_dim, is_log_integrand, 1);

        double * tuning_sample_points = new double[tuning_num_points * truncated_ndim];

        slice_sampler.Initialize();

        slice_sampler.sample(&posterior, posterior.mode_point,
                             initial_sampling_range, 1, tuning_sample_points,
                             tuning_num_points);
        
        best_autocor_offset =
            best_autocorrelation_offset(tuning_sample_points, truncated_ndim, 
                                        tuning_num_points, autocor_max_offset, 
                                        autocor_max_value);
        if (verbose)
        {
            fprintf(stdout, "SS: %Zu, position: %Zu\n", 
                    best_autocor_offset, locus.position);
        }
        
        delete tuning_sample_points;


        if (best_autocor_offset <= target_autocor_offset)
        {
            //slice sampling succeeded
            double * slice_sample_points = new double[final_num_points * truncated_ndim];
            double * sample_points_buf = new double[final_num_points * 4];
            std::vector<double *> sample_points_sortable(final_num_points);

            slice_sampler.sample(&posterior, posterior.mode_point,
                                 initial_sampling_range, 
                                 best_autocor_offset,
                                 slice_sample_points,
                                 final_num_points);

            add_normalized_dimension(slice_sample_points, truncated_ndim,
                                     final_num_points, sample_points_buf,
                                     &sample_points_sortable);
            
            sprintf(line_label, "%s\t%s\t%s\t%Zu\t%c\t%Zu\t%Zu\t%5.5lf\t%5.5lf", label_string, 
                    "SS", locus.reference, 
                    locus.position, 
                    locus.reference_base, 
                    locus.read_depth, effective_depth,
                    pos_anomaly_score, neg_anomaly_score);
            
            //approximate marginal modes (expensive!)
//             std::fill(marginal_modes, marginal_modes + 4, -1.0);
//             if (marginal_mode_num_points > 0)
//             {
//                 for (size_t d = 0; d != full_ndim; ++d)
//                 {
//                     marginal_modes[d] = 
//                         window_averaged_mode(& sample_points_sortable, d, 
//                                              marginal_mode_num_points);
//                 }
//             }
            
            print_quantiles(posterior_output_fh, & sample_points_sortable, 
                            posterior.mode_point, line_label, 
                            dimension_labels, "+", quantiles, num_quantiles,
                            truncated_ndim + 1);


            if (cdfs_output_fh != NULL)
            {
                print_numerical_cdfs(cdfs_output_fh, line_label, & sample_points_sortable, full_ndim);
            }

            fflush(posterior_output_fh);
            delete sample_points_buf;
            delete slice_sample_points;

            continue;
        }


        //nothing else worked.
    }
    
    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
    {
        fclose(cdfs_output_fh);
    }

    gsl_rng_free(rand_gen);

    delete prior_alphas;
    delete quantiles;
    delete data_reader;
    return 0;

}
