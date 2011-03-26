#include <cstdio>
#include <cstdlib>

int main(int argc, char ** argv){
    
    FILE * first;
    FILE * second;

    if (argc != 3)
    {
        fprintf(stderr, "Usage: shuffle_paired_fastq read1.fq read2.fq > reads_shuffled.fq");
        return 1;
    }

    first = fopen(argv[1], "r");
    second = fopen(argv[2], "r");
    
    if (first == NULL)
    {
        fprintf(stderr, "Couldn't open first fastq file %s\n",
                argv[1]);
        exit(2);
    }

    if (second == NULL)
    {
        fprintf(stderr, "Couldn't open second fastq file %s\n",
                argv[2]);
        exit(2);
    }

    char id[1024];
    char sequence[4096];
    char spacer[1024];
    char quality_string[4096];

    while (fscanf(first, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string) == 4)
    {
        fprintf(stdout, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);

        fscanf(second, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);
        fprintf(stdout, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);

    }
    fclose(first);
    fclose(second);
    return 0;
}
