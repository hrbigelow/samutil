#include "sim_expression.h"
#include "transcript_generator.h"

#include <cstdio>
#include <cstring>

int sim_expression_usage()
{
    fprintf(stderr,
            "Usage:\nsim expression density_function_string transcripts.gtf transcript_expression.txt\n\n"
            "transcript_expression.txt output format:\ntranscript_id\tlength\tsimulated_expression\n\n"
            "Supported density function distributions:\n\n%s\n",
            supported_distributions);
    return 1;
}

int main_sim_expression(int argc, char **argv)
{

    char c;
    while ((c = getopt(argc, argv, "")) >= 0)
    {
        switch(c)
        {
        default: return sim_expression_usage(); break;
        }
    }

    int req_args = 3;
    int req_arg_count = optind + req_args;

    if (argc != req_arg_count)
    {
        return sim_expression_usage();
    }

    char * density_function = argv[optind];
    char * transcript_gtf_file = argv[optind + 1];
    char * transcript_expression_file = argv[optind + 2];

    TranscriptGenerator transcript_gen;
    transcript_gen.Initialize(transcript_gtf_file, "exon");

    transcript_gen.GenerateExpression(density_function);

    transcript_gen.PrintExpression(transcript_expression_file);

    return 0;
}
