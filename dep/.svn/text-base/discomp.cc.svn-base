#include "tools.h"
#include "error_estimate.h"
#include "integrands.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "simulation.h"
#include "nucleotide_reader.h"
#include "nucleotide_stats.h"
#include "base_qual_strand_reader.h"

#include "usage_strings.h"

int discomp_usage()
{
    fprintf(stderr, 
            "Usage: dep discomp [options] jpd_paramss input.pileup output.comp\n" 
            "Options:\n\n"
            "-d STRING   %s\n"
            "-l STRING   %s\n"
            "-e FILE     file containing list of compositions to evaluate posterior\n"
            "-q INT      %s\n"
            "\n",
            Usage::data_adaptor, Usage::label_string, Usage::quality_string);
    return 1;
}


int main_discomp(int argc, char ** argv)
{

    char input_data_type[1000];
    strcpy(input_data_type, "base_qual_strand");
    
    char label_string[1000];
    char evaluation_points_file[1000];
    size_t min_quality_score = 5;

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* posterior_output_file;

    char c;
    while ((c = getopt(argc, argv, "d:l:e:q:")) >= 0)
    {
        switch(c)
        {
        case 'd': strcpy(input_data_type, optarg); break;
        case 'l': strcpy(label_string, optarg); break;
        case 'e': strcpy(evaluation_points_file, optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        }
    }
    if (argc - optind != 3)
    {
        return discomp_usage();
    }

    jpd_data_params_file = argv[optind + 1];
    pileup_input_file = argv[optind + 2];
    posterior_output_file = argv[optind + 3];

    size_t num_evaluation_points;

    printf ("input_data_type: %s\n"
            "label_string: %s\n"
            "evaluation_points_file: %s\n"
            "min_quality_score: %Zu\n"
            "jpd_data_params_file: %s\n"
            "pileup_input_file: %s\n"
            "posterior_output_file: %s\n",

            input_data_type,
            label_string,
            evaluation_points_file,
            min_quality_score,
            jpd_data_params_file,
            pileup_input_file,
            posterior_output_file
            );
    

    double * evaluation_points;
    double * evaluation_point_dist;

    double * numbers_buffer = 
        parse_basecomp_prior_file(evaluation_points_file,
                                  &evaluation_points, 
                                  &evaluation_point_dist);

    
    FILE * posterior_output_fh = fopen(posterior_output_file, "w");
    if (posterior_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open posterior_output_file %s\n",
                posterior_output_file);
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
    NucleotideStats global_params = data_reader->read_from_rdb(jpd_data_params_file);
    
    //we are integrating the actual posterior
    //char const* dimension_labels[] = { "A", "C", "G", "T" };
    //char line_label[1000];

    size_t effective_depth;

    //double composition_prior_alphas[] = { 1e-5, 1e-5, 1e-5, 1e-5 };

    ErrorEstimate default_model;
    default_model.set_model_params(& global_params);
    default_model.set_discrete_prior_dist(evaluation_points, 
                                          evaluation_point_dist, 
                                          num_evaluation_points);

    double * posterior_values = new double[num_evaluation_points];
    double * log_posterior_values = new double[num_evaluation_points];

    while (data_reader->more_loci())
    {
        LocusSummary locus = 
            data_reader->get_next_locus(global_params, 
                                        static_cast<void const*>(& min_quality_score));

        //set default model posterior
        default_model.set_locus_data(& locus);

        effective_depth = locus.read_depth;

        //evaluate posterior at all points consistent with ploidy
        double * evaluation_point;
        for (size_t p = 0; p != num_evaluation_points; ++p)
        {
            evaluation_point = evaluation_points + (p * 4);
            log_posterior_values[p] = 
                default_model.log_likelihood(evaluation_point)
                + default_model.log_dirichlet_prior(evaluation_point);
        }

        double * max_val = std::max_element(log_posterior_values, 
                                            log_posterior_values + num_evaluation_points);
        size_t max_elem = std::distance(log_posterior_values, max_val);

        for (size_t p = 0; p != num_evaluation_points; ++p)
        {
            posterior_values[p] = exp(log_posterior_values[p] - *max_val);
        }

        fprintf(posterior_output_fh, "%s\t%s\t%Zu\t%c\t%Zu\t%Zu\t%Zu", 
                label_string, locus.reference, 
                locus.position, locus.reference_base, 
                locus.read_depth, effective_depth,
                max_elem);
        
        normalize(posterior_values, num_evaluation_points, posterior_values);

        //adjust values, normalize
        for (size_t p = 0; p != num_evaluation_points; ++p)
        {
            fprintf(posterior_output_fh, "\t%g", posterior_values[p]);
        }

        fprintf(posterior_output_fh, "\n");

    }
    fclose(posterior_output_fh);

    delete data_reader;
    delete numbers_buffer;

    delete posterior_values;
    delete log_posterior_values;
    return 0;

}
