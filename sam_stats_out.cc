#include "sam_stats.h"
#include "sam_stats_aux.h"
#include "matrix_tools.h"

#include <map>
#include <cstdio>

//process the raw stats output file and produce human-usable statistics
int main_out(int argc, char ** argv)
{

    if (argc != 2)
    {
        return out_usage();
    }

    char * raw_stats_file = argv[1];

    size_t num_consensus_match;
    size_t num_strands;
    size_t num_consensus_scores;
    size_t num_quality_scores;
    size_t num_alignment_scores;
    size_t num_read_positions;

    FILE * raw_stats_fh = fopen(raw_stats_file, "r");
    if (raw_stats_fh == NULL)
    {
        fprintf(stderr, "Couldn't open raw statistics file %s\n", raw_stats_file);
        exit(1);
    }

    fscanf(raw_stats_fh, "%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n",
           &num_read_positions,
           &num_consensus_match,
           &num_consensus_scores,
           &num_quality_scores,
           &num_strands,
           &num_alignment_scores);

    std::map<int, int> alignment_scores;
    std::map<int, int>::const_iterator score_iter;

    int score;
    int index;
    for (size_t a = 0; a != num_alignment_scores; ++a)
    {
        fscanf(raw_stats_fh, "%i\t%i\n", &index, &score);
        alignment_scores[score] = index;
    }

    size_t const ndims = 6;

    //dimensions (lowest first): RP CM CS QS ST AS
    size_t dim_sizes[ndims];
    dim_sizes[0] = num_read_positions;
    dim_sizes[1] = num_consensus_match;
    dim_sizes[2] = num_consensus_scores;
    dim_sizes[3] = num_quality_scores;
    dim_sizes[4] = num_strands;
    dim_sizes[5] = num_alignment_scores;


    size_t num_categories = 1;
    for (size_t d = 0; d != ndims; ++d)
    {
        num_categories *= dim_sizes[d];
    }

    int * alignment_stats = new int[num_categories];
    fread(alignment_stats, sizeof(int), num_categories, raw_stats_fh);
    fclose(raw_stats_fh);

    //Compile and print all statistics

    // all single marginals
    char const* tags[] = {
        "read_positions",
        "consensus_match",
        "consensus_scores",
        "quality_scores",
        "strands",
        "alignment_scores"
    };

    for (size_t d = 0; d != 6; ++d)
    {
        int * marginal_dist = new int[dim_sizes[d]];
        std::fill(marginal_dist, marginal_dist + dim_sizes[d], 0);
        
        //tally alignment_score_dist
        bool marginal_mask[] = { false, false, false, false, false, false };

        marginal_mask[d] = true;

        MatrixTools::marginalize(alignment_stats, dim_sizes, ndims, 
                                 marginal_mask, marginal_dist);
        
        for (size_t i = 0; i != dim_sizes[d]; ++i)
        {
            if (marginal_dist[i] > 0)
            {
                fprintf(stdout, "%s\t%Zu\t%i\n", tags[d], i, marginal_dist[i]);
            }
        }
        delete marginal_dist;
    }
    
    //1. AS

    /*
    int * alignment_score_dist = new int[num_alignment_scores];
    std::fill(alignment_score_dist, alignment_score_dist + num_alignment_scores, 0);

    //tally alignment_score_dist
    bool alignment_score_marg_mask[] = { false, false, false, false, false, true };

    MatrixTools::marginalize(alignment_stats, dim_sizes, ndims, 
                             alignment_score_marg_mask, alignment_score_dist);
    
    for (score_iter = alignment_scores.begin();
         score_iter != alignment_scores.end(); ++score_iter)
    {
        int score = (*score_iter).first;
        int index = (*score_iter).second;
        if (alignment_score_dist[index] > 0)
        {
            fprintf(stdout, "alignment_score\t%i\t%i\n", score, alignment_score_dist[index]);
        }
    }
    delete alignment_score_dist;
    */

    //2. consensus_match X consensus_scores X quality_scores X alignment_scores
    //2. CM x CS x QS x AS
    bool calibration_marg_mask[] = { false, true, true, true, false, true };

    size_t num_calibration = 
        dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[5];

    int * calibration_dist = new int[num_calibration];

    MatrixTools::marginalize(alignment_stats, dim_sizes, ndims,
                             calibration_marg_mask, calibration_dist);
    

    size_t calibration_cat[] = { 0, 0, 0, 0 };
    size_t calibration_dim_sizes[] = {
        dim_sizes[1],
        dim_sizes[2],
        dim_sizes[3],
        dim_sizes[5]
    };

    int movable_digits[] = { 0, 1, 2, 3, -1 };

    for (size_t merged_index = 0; merged_index != num_calibration; ++merged_index)
    {
        if (calibration_dist[merged_index] > 0)
        {
            fprintf(stdout, "CM_X_CS_X_QS_X_AS\t%Zu\t%Zu\t%Zu\t%i\t%i\n", 
                    calibration_cat[0],
                    calibration_cat[1],
                    calibration_cat[2],
                    alignment_scores[calibration_cat[3]],
                    calibration_dist[merged_index]);
        }
        MatrixTools::increment_counter(calibration_dim_sizes, movable_digits,
                                       calibration_cat);
    }
                    
    delete calibration_dist;

    //3. What else?

    
    delete alignment_stats;

    return 0;
}
