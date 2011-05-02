#ifndef _TRANSCRIPT_GENERATOR_H
#define _TRANSCRIPT_GENERATOR_H

#include <cstdlib>
#include <string>
#include <set>

#include <gsl/gsl_rng.h>

#include "gtf.h"

class TranscriptGenerator
{
    gsl_rng * rand_gen;

public:
    TRANSCRIPT_EXPR transcript_info;
    TRANSCRIPT_EXPR::iterator tr_iter;

    //initialize various model parameters for generating the transcript levels
    TranscriptGenerator();
    ~TranscriptGenerator();

    void Initialize(char const* transcript_file,
                    char const* feature);

    void GenerateExpression(char const* density_function_string, size_t random_seed);
    void PrintExpression(char const* expression_file);

};

extern char const* supported_distributions;

#endif // _TRANSCRIPT_GENERATOR_H
