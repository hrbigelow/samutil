#include <vector>
#include <utility>
#include <numeric>

#include <gsl/gsl_math.h>

#include "tools.h"
#include "error_estimate.h"
#include "integrands.h"
#include "metropolis.h"
#include "sampling.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "stats_tools.h"
#include "nucleotide_reader.h"
#include "nucleotide_stats.h"
#include "base_qual_strand_reader.h"

#include "usage_strings.h"

int anomaly_usage()
{
    fprintf(stderr, 
            "\nUsage: dep anomaly [options] jpd_params input.pileup output.anom\n" 
            "Options:\n\n"
            "-d STRING   %s\n"
            "-l STRING   %s\n"
            "-q INT      %s\n"
            "\n", Usage::data_adaptor, Usage::label_string, Usage::quality_string
            );
    return 1;
}


int main_anomaly(int argc, char ** argv)
{

    char input_data_type[100];
    strcpy(input_data_type, "base_qual_strand");

    char label_string[100];
    strcpy(label_string, "comp");

    int min_quality_score = 5;

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* anomaly_output_file;

    char c;
    while ((c = getopt(argc, argv, "d:l:q:")) >= 0)
    {
        switch(c)
        {
        case 'd': strcpy(input_data_type, optarg); break;
        case 'l': strcpy(label_string, optarg); break;
        case 'q': min_quality_score = atoi(optarg); break;
        }
    }
    if (argc - optind < 3)
    {
        return anomaly_usage();
    }

    jpd_data_params_file = argv[optind + 1];
    pileup_input_file = argv[optind + 2];
    anomaly_output_file = argv[optind + 3];

    
    size_t full_ndim = 4;
    size_t truncated_ndim = 3;

    FILE * anomaly_output_fh = fopen(anomaly_output_file, "w");
    if (anomaly_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open anomaly_output_file %s\n",
                anomaly_output_file);
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
    char const* dimension_labels[] = { "A", "C", "G", "T" };
    char line_label[1000];

    double mode_tolerance = 1e-60;
    size_t max_modefinding_iterations = 3000;

    bool use_independence_chain_mh = true;
    bool may_underflow = true;

    size_t effective_depth;
    
    ErrorEstimate default_model;
    default_model.set_model_params(& global_params);

    while (data_reader->more_loci())
    {
    }

    fclose(anomaly_output_fh);
    delete data_reader;
    return 0;

}
