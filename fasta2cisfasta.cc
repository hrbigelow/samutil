#include <cstdio>
#include <cstdlib>

#include "dep/tools.h"

int main(int argc, char ** argv){

    char const* input_fasta_file = argv[1];
    char const* output_fasta_file = argv[2];

    if (argc != 3)
    {
        fprintf(stderr, "removes newlines from fasta file\n\n");
        fprintf(stderr, "usage: fasta2cisfasta input.fa output.fa\n\n");
        return 1;
    }

    FILE * input_fasta_fh = open_or_die(input_fasta_file, "r", "Input fasta file");
    FILE * output_fasta_fh = open_or_die(output_fasta_file, "w", "Converted output fasta file");

    const int MAX_LINE=100000;
    int chunk = 0;
    char * line = new char[MAX_LINE];
    
    while (! feof(input_fasta_fh)){
        
        fscanf(input_fasta_fh, "%s", line);
        if (line[0] == '>'){ // we have a header line
            ++chunk;
            if (chunk > 1)
            { 
                fprintf(output_fasta_fh, "\n");
            }
            fprintf(output_fasta_fh, "%s\n", line);
        } 
        else
        {  // print the DNA chunk without a newline
            fprintf(output_fasta_fh, "%s", line); 
        }
    }
    
    fprintf(output_fasta_fh, "\n");

    fclose(input_fasta_fh);
    fclose(output_fasta_fh);

    delete line;
    
    return 0;
}

