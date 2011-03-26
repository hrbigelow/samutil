#ifndef _FRAGMENT_GENERATOR_H
#define _FRAGMENT_GENERATOR_H

#include <vector>
#include <utility>
#include <map>
#include <string>

#include <gsl/gsl_rng.h>

#include "gtf.h"



class FragmentGenerator
{

 private:
    gsl_rng * rand_gen;

 public:
    enum { CUT, UNIFORM } sampling_scheme;
    TRANSCRIPT_EXPR transcript_info;

    //for 'UNIFORM' mode
    size_t num_fragment_lengths;
    size_t * fragment_lengths;
    size_t step_size;

    //for 'CUT' mode
    size_t fragment_length_min;
    size_t fragment_length_max;
    float cut_probability;
    size_t total_transcript_molecules;

    FILE * transcript_expression_fh;

 public:
    typedef std::vector<std::pair<size_t, size_t> > BOUNDS;


    
    FragmentGenerator();
    ~FragmentGenerator();

    void Initialize(char const* scheme_string,
                    char const* expression_file);

    BOUNDS Sample(unique_transcript const& transcript);


    //generate uniform coverage across a transcript, for a given step
    //size, and a series of fragment lengths in [frag_length_start,
    //frag_length_end)
    BOUNDS Uniform(size_t transcript_length);
    
    //generates a series of random cuts in the transcript molecule, storing
    //only the fragments of acceptable length
    BOUNDS Cut(size_t transcript_length, size_t num_transcript_mols);
};

#endif // _FRAGMENT_GENERATOR_H
