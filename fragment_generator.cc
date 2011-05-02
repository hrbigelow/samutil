#include "fragment_generator.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cstring>
#include <cassert>
#include <algorithm>

FragmentGenerator::FragmentGenerator()
{
    this->rand_gen = gsl_rng_alloc(gsl_rng_taus);
}


void FragmentGenerator::Initialize(char const* scheme_string,
                                   char const* expression_file,
                                   size_t random_seed)
{
    gsl_rng_set(this->rand_gen, random_seed);

    FILE * expression_fh = fopen(expression_file, "r");
    if (expression_fh == NULL)
    {
        fprintf(stderr, "Error: FragmentGenerator::Initialize(): "
                "Couldn't open transcript expression file %s\n",
                expression_file);
        exit(1);
    }
    //parse the transcript expression file, format is:
    //transcript_id length expression_level
    float expression_level;
    size_t transcript_length;
    char transcript_id[1024];
    char contig_name[1024];
    char strand;

    while (! feof(expression_fh))
    {
        fscanf(expression_fh, "%s\t%s\t%c\t%zu\t%f\n", transcript_id, contig_name, &strand, 
               &transcript_length, &expression_level);

        this->transcript_info[unique_transcript(contig_name, transcript_id, strand)] = 
            TranscriptInfo(expression_level, 0, transcript_length);
    }
    fclose(expression_fh);

    size_t fragment_length_step;

    if (strncmp(scheme_string, "cut", 3) == 0)
    {
        this->sampling_scheme = CUT;
        int parse = sscanf(scheme_string, "cut:%f,%zu,%zu,%zu",
                           &this->cut_probability, 
                           &this->fragment_length_min,
                           &this->fragment_length_max,
                           &this->total_transcript_molecules);
        if (parse != 4)
        {
            fprintf(stderr, "Error: FragmentGenerator::Initialize(): "
                    "invalid format for cut sampling\n"
                    "Must be 'cut:%%f,%%i,%%i,%%i'\n");
            exit(1);
        }

        if (! (this->fragment_length_min <= this->fragment_length_max))
        {
            fprintf(stderr, "FragmentGenerator::Initialize(): "
                    "provided min fragment length (%Zu) is not less than or equal to "
                    "provided max fragment length (%Zu)\n",
                    this->fragment_length_min, this->fragment_length_max);
            exit(1);
        }

        fragment_length_step = 1;
        this->num_fragment_lengths = this->fragment_length_max - this->fragment_length_min + 1;
        this->step_size = 0;

        //generate transcript molecules
        size_t multinomial_K = this->transcript_info.size();

        TRANSCRIPT_EXPR::iterator ti_iter;

        double * multinomial_probs = new double[multinomial_K];
        size_t t;
        for (t = 0, ti_iter = this->transcript_info.begin(); 
             t != multinomial_K; ++t, ++ti_iter)
        {
            multinomial_probs[t] = (*ti_iter).second.relative_level;
        }


        unsigned int * transcript_counts = new unsigned int[multinomial_K];

        gsl_ran_multinomial(this->rand_gen, multinomial_K, 
                            this->total_transcript_molecules, 
                            multinomial_probs, transcript_counts);
        
        for (t = 0, ti_iter = this->transcript_info.begin(); 
             t != multinomial_K; ++t, ++ti_iter)
        {
            (*ti_iter).second.count = transcript_counts[t];
        }
        delete multinomial_probs;
        delete transcript_counts;
    }

    else if (strncmp(scheme_string, "uniform", 7) == 0)
    {
        this->sampling_scheme = UNIFORM;
        int parse = sscanf(scheme_string, "uniform:%zu,%zu,%zu,%zu", 
                           &this->step_size, 
                           &this->fragment_length_min, 
                           &fragment_length_step,
                           &this->fragment_length_max);

        if (parse != 4)
        {
            fprintf(stderr, "Error: FragmentGenerator::Initialize(): "
                    "invalid format for uniform sampling\n"
                    "Must be 'uniform:%%i,%%i,%%i,%%i'\n");
            exit(1);
        }

        if (this->fragment_length_min == this->fragment_length_max)
        {
            this->num_fragment_lengths = 1;
        }
        else
        {
            this->num_fragment_lengths = 1 + 
                (this->fragment_length_max - this->fragment_length_min) / fragment_length_step;

            assert((this->fragment_length_max - this->fragment_length_min) % fragment_length_step == 0);
        }

    }
    else
    {
        fprintf(stderr, "Error: unknown sampling scheme, %s\n"
                "Should be one of 'cut:' or 'uniform:'\n", scheme_string);
        exit(1);
    }

    size_t len, fi;
    this->fragment_lengths = new size_t[this->num_fragment_lengths];
    for (len = this->fragment_length_min, fi = 0; 
         fi != this->num_fragment_lengths; len += fragment_length_step, ++fi)
    {
        this->fragment_lengths[fi] = len;
    }
    
}


FragmentGenerator::~FragmentGenerator()
{
    gsl_rng_free(this->rand_gen);
    delete this->fragment_lengths;
}


FragmentGenerator::BOUNDS 
FragmentGenerator::Sample(unique_transcript const& transcript)
{
    TRANSCRIPT_EXPR::iterator tit;
    tit = this->transcript_info.find(transcript);
    if (tit == this->transcript_info.end())
    {
        fprintf(stderr, "Error: FragmentGenerator::Sample: couldn't find"
                " transcript information for transcript %s on %s (%c)\n",
                transcript.transcript_id.c_str(), 
                transcript.contig_name.c_str(), 
                transcript.strand);

        exit(1);
    }
    TranscriptInfo & tinfo = (*tit).second;

    switch(this->sampling_scheme)
    {
    case UNIFORM: return this->Uniform(tinfo.length); break;
    case CUT: return this->Cut(tinfo.length, tinfo.count); break;
    }
}

FragmentGenerator::BOUNDS FragmentGenerator::Uniform(size_t transcript_length)
{
    FragmentGenerator::BOUNDS bounds;
    
    for (size_t f = 0; f != this->num_fragment_lengths; ++f)
    {
        if (transcript_length < this->fragment_lengths[f])
        {
            continue;
        }
        size_t max_start_pos = transcript_length - this->fragment_lengths[f];
        for (size_t start_pos = 0; start_pos <= max_start_pos; start_pos += this->step_size)
        {
            bounds.push_back(std::make_pair(start_pos, start_pos + this->fragment_lengths[f]));
        }
    }
    return bounds;
}


FragmentGenerator::BOUNDS 
FragmentGenerator::Cut(size_t transcript_length, size_t num_transcript_mols)
{
    FragmentGenerator::BOUNDS bounds;
    for (size_t t = 0; t != num_transcript_mols; ++t)
    {
        size_t num_cuts = 
            gsl_ran_binomial(this->rand_gen, this->cut_probability, 
                             transcript_length - 1) + 2;

        size_t * cut_positions = new size_t[num_cuts];
        cut_positions[0] = 0;
        cut_positions[num_cuts - 1] = transcript_length;

        for (size_t cut = 1; cut != num_cuts - 1; ++cut)
        {

            //gsl_rng_uniform_int returns random int in range [0,n-1] inclusive
            cut_positions[cut] = 
                gsl_rng_uniform_int(this->rand_gen, transcript_length - 1 - cut) + 1;

            if (cut > 1)
            {
                for (size_t prev_cut = 0; prev_cut != cut; ++prev_cut)
                {
                    if (cut_positions[cut] >= cut_positions[prev_cut])
                    {
                        ++cut_positions[cut];
                    }
                }
            }
        }
        std::sort(cut_positions, cut_positions + num_cuts);
        for (size_t frag = 1; frag != num_cuts; ++frag)
        {
            size_t frag_length = cut_positions[frag] - cut_positions[frag - 1];
            if (this->fragment_length_min <= frag_length
                && frag_length <= this->fragment_length_max)
            {
                bounds.push_back(std::make_pair(cut_positions[frag - 1], cut_positions[frag]));
            }
        }
        delete cut_positions;
    }
    return bounds;
}
