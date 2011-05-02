#include "transcript_generator.h"
#include "gtf.h"

#include "dep/tools.h"

#include <cstdio>
#include <cstring>

#include <gsl/gsl_randist.h>

TranscriptGenerator::TranscriptGenerator()
{
    this->rand_gen = gsl_rng_alloc(gsl_rng_taus);
}

TranscriptGenerator::~TranscriptGenerator()
{
    gsl_rng_free(this->rand_gen);
}


//tally exon length from each transcript
void TranscriptGenerator::Initialize(char const* gtf_file, 
                                     char const* feature)
{
    FILE * gtf_fh = fopen(gtf_file, "r");
    if (gtf_fh == NULL)
    {
        fprintf(stderr, "Couldn't open gtf file %s\n", gtf_file);
        exit(1);
    }

    GTFEntry gtf_entry;

    std::map<unique_transcript, size_t> transcript_lengths;
    std::map<unique_transcript, size_t>::iterator tliter;

    unique_transcript transcript;

    while (gtf_entry.get_next_record(gtf_fh))
    {
        if (strcmp(feature, gtf_entry.feature) != 0)
        {
            continue;
        }
        transcript = unique_transcript(gtf_entry.seqname, gtf_entry.transcript_id,
                                       gtf_entry.strand);

        tliter = transcript_lengths.find(transcript);
        if (tliter == transcript_lengths.end())
        {
            tliter = transcript_lengths.insert(transcript_lengths.end(),
                                               std::make_pair(transcript, 0));
        }
        (*tliter).second += gtf_entry.length();
    }

    for (tliter = transcript_lengths.begin(); tliter != transcript_lengths.end(); ++tliter)
    {
        this->transcript_info.insert(std::make_pair
                                     ((*tliter).first,
                                      TranscriptInfo(0.0, 0, (*tliter).second)));
    }
    fclose(gtf_fh);
}


//generates expression using a density function string to describe the function
//and parameters
char const* supported_distributions = 
    "lognormal:%%f,%%f                        zeta, sigma\n"
    "tln:%%f,%%f,%%f   (truncated log-normal) fraction_expressed, zeta, sigma\n";


void TranscriptGenerator::GenerateExpression(char const* density_function_string,
                                             size_t random_seed)
{
    gsl_rng_set(this->rand_gen, random_seed);

    std::set<TranscriptInfo const*>::const_iterator it;

    double zeta;
    double sigma;
    double fraction_expressed;
    if (sscanf(density_function_string, "lognormal:%lf,%lf", &zeta, &sigma) == 2)
    {
        for (TRANSCRIPT_EXPR::iterator it = this->transcript_info.begin(); 
             it != this->transcript_info.end(); ++it)
        {
            (*it).second.relative_level = gsl_ran_lognormal(this->rand_gen, zeta, sigma);
        }
    }
    else if (sscanf(density_function_string, "tln:%lf,%lf,%lf", 
                    &fraction_expressed, &zeta, &sigma) == 3)
    {
        double expressed;
        for (TRANSCRIPT_EXPR::iterator it = this->transcript_info.begin(); 
             it != this->transcript_info.end(); ++it)
        {
            expressed = gsl_rng_uniform(this->rand_gen);
            
            (*it).second.relative_level = (expressed < fraction_expressed) 
                ? gsl_ran_lognormal(this->rand_gen, zeta, sigma)
                : 0L;
        }
        
    }
    else
    {
        fprintf(stderr, "TranscriptGenerator::GenerateExpression: "
                "Error: don't understand this density function description: '%s'\n\n"
                "Must be one of the following:\n\n%s\n",
                density_function_string, supported_distributions);
        exit(1);
    }
}


void TranscriptGenerator::PrintExpression(char const* expression_file)
{
    FILE * expression_fh = open_if_present(expression_file, "w");
    TRANSCRIPT_EXPR::iterator it;
    for (it = this->transcript_info.begin(); it != this->transcript_info.end(); ++it)
    {
        unique_transcript const& ut = (*it).first;
        TranscriptInfo const& ti = (*it).second;

        fprintf(expression_fh, "%s\t%s\t%c\t%Zu\t%f\n", 
                ut.transcript_id.c_str(), ut.contig_name.c_str(), ut.strand,
                ti.length, ti.relative_level);
    }
    fclose(expression_fh);
}
