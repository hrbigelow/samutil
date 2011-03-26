#include "align_eval_stats.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>

#include "dep/tools.h"

int align_eval_stats_usage()
{
    fprintf(stderr,
            "Usage:\n\n"
            "align_eval stats [OPTIONS] jumps.txt stats.txt dist.txt\n\n"
            "Options:\n\n");

    fprintf(stderr,
            "jumps.txt has fields:\n"
            "<contig> <position> <guide_jump> <correct_jump> <error_jump>\n"
            "See 'align_eval raw' to produce jumps.txt file from SAM alignment\n\n"
            "dist.txt has fields <fraction> <num_correct> <num_correct_cumul> <num_error> <num_error_cumul>\n"
            "<fraction> is the fraction of correct or error\n"
            "stats.txt has fields <guide_depth> <correct_depth> <error_depth> <locus_count>\n\n");
    return 1;
}


struct Depths
{
    int guide : 20;
    int correct : 20;
    int error : 24;
};

union TripleDepth
{
    uint64_t index;
    Depths depth;

    bool operator<(TripleDepth const& b) const
    {
        return this->index < b.index;
    }
};


struct margin_cumul
{
    size_t correct;
    size_t correct_cumul;
    size_t error;
    size_t error_cumul;
};
        



int main_align_eval_stats(int argc, char ** argv)
{
    char c;
    while ((c = getopt(argc, argv, "")) >= 0)
    {
        switch(c)
        {
        default: return align_eval_stats_usage(); break;
        }
    }
    
    int arg_count = optind + 3;
    if (argc != arg_count)
    {
        return align_eval_stats_usage();
    }

    char const* jumps_file = argv[optind];
    char const* stats_file = argv[optind + 1];
    char const* dist_file = argv[optind + 2];

    FILE * jumps_fh = open_or_die(jumps_file, "r", "Input jumps file");
    FILE * stats_fh = open_or_die(stats_file, "w", "Output stats file");
    FILE * dist_fh = open_or_die(dist_file, "w", "Output dist file");

    //maintain a running current total of guide, correct, and error.
    char contig[256] = "";
    char prev_contig[256] = "";
    size_t bound = 0;
    size_t prev_bound;
    int guide_jump, correct_jump, error_jump;

    std::map<TripleDepth, size_t> depth_count;
    std::map<TripleDepth, size_t>::const_iterator dit;

    TripleDepth triple;
    triple.depth.guide = 0;
    triple.depth.correct = 0;
    triple.depth.error = 0;

    /*
      idea: maintain a set of current depths (guide, correct, error).
      at each jump, record the number of loci between the last jump and this one,
      and increment the count of current depths for 
     */

    while (! feof(jumps_fh))
    {
        //strcpy(prev_contig, contig);

        fscanf(jumps_fh, "%s\t%zu\t%i\t%i\t%i\n",
               contig, &bound, &guide_jump, &correct_jump, &error_jump);

        if (strcmp(prev_contig, contig) != 0)
        {
            //on a new contig.  
            assert(triple.depth.guide == 0);
            assert(triple.depth.correct == 0);
            assert(triple.depth.error == 0);
            prev_bound = 0;
            //what to do about knowing the end of the contig?
        }
        dit = depth_count.find(triple);
        if (dit == depth_count.end())
        {
            depth_count[triple] = 0;
        }
        depth_count[triple] += (bound - prev_bound);

        triple.depth.guide += guide_jump;
        triple.depth.correct += correct_jump;
        triple.depth.error += error_jump;

        prev_bound = bound;
        strcpy(prev_contig, contig);
    }

    assert(triple.depth.guide == 0);
    assert(triple.depth.correct == 0);
    assert(triple.depth.error == 0);

    //compute distributions marginalized over % correct and % error
    fprintf(stats_fh, "#%s\t%s\t%s\t%s\n", 
            "num_simulated_bases", 
            "num_correct_bases",
            "num_error_bases", 
            "num_loci");

    std::map<float, margin_cumul> margins;

    float fraction;
    size_t total_num_loci = 0;
    for (dit = depth_count.begin(); dit != depth_count.end(); ++dit)
    {
        TripleDepth const& t = (*dit).first;
        size_t c = (*dit).second;
        fprintf(stats_fh, "%i\t%i\t%i\t%Zu\n", t.depth.guide, t.depth.correct, t.depth.error, c);
        if (t.depth.guide > 0)
        {
            fraction = static_cast<float>(t.depth.correct) / static_cast<float>(t.depth.guide);
            margins[fraction].correct += c;

            fraction = static_cast<float>(t.depth.error) / static_cast<float>(t.depth.guide);
            margins[fraction].error += c;
            total_num_loci += c;
        }
    }

    std::map<float, margin_cumul>::iterator mit;
    size_t correct_cumul = 0;
    size_t error_cumul = total_num_loci;

    fprintf(dist_fh, "#%s\t%s\t%s\t%s\t%s\n", "fraction_of_guide_depth", "num_correct_loci",
            "num_correct_loci_better", "num_error_loci", "num_error_loci_better");

    for (mit = margins.begin(); mit != margins.end(); ++mit)
    {
        margin_cumul & d = (*mit).second;
        d.correct_cumul = correct_cumul;
        d.error_cumul = error_cumul;
        correct_cumul += d.correct;
        error_cumul -= d.error;
        fprintf(dist_fh, "%f\t%Zu\t%Zu\t%Zu\t%Zu\n", (*mit).first,
                d.correct, d.correct_cumul, d.error, d.error_cumul);
    }

    assert(error_cumul == 0);

    fclose(stats_fh);
    fclose(dist_fh);
    fclose(jumps_fh);

    return 0;
}
