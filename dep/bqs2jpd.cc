#include "nucleotide_stats.h"
#include "tools.h"


int bqs2jpd_usage()
{
    fprintf(stderr,
            "\nUsage: dep bqs2jpd counts.bqs > jpd_priors\n");
    return 1;
}

int main_bqs2jpd(int argc, char ** argv)
{

    if (argc != 2)
    {
        return bqs2jpd_usage();
    }

    char basecall;
    size_t quality;
    char strand;
    float count;
    float error_prob;
    float data_prob[4];
    size_t basecall_index;

    FILE * counts_fh = fopen(argv[1], "r");

    while (! feof(counts_fh))
    {
        fscanf(counts_fh, "%c\t%zi\t%c\t%f\n", &basecall, &quality, &strand, &count);
        error_prob = QualityToErrorProb(quality);

        std::fill(data_prob, data_prob + 4, (error_prob / 3.0) * count);

        basecall_index = 
            Nucleotide::base_to_index[static_cast<size_t>(basecall)];

        data_prob[basecall_index] = (1.0 - error_prob) * count;
        
        fprintf(stdout, "%c_%Zu_%c\t%8.2f\t%8.2f\t%8.2f\t%8.2f\n",
                basecall, quality, strand, data_prob[0], data_prob[1],
                data_prob[2], data_prob[3]);
    }

    fclose(counts_fh);

    return 0;
}
