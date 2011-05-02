#include <cstdio>
#include <cstdlib>

#include "dep/tools.h"
#include "file_utils.h"

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
    char * line = new char[MAX_LINE + 1];
    char next_ch;
    while (! feof(input_fasta_fh)){
        
        next_ch = fgetc(input_fasta_fh);
        ungetc(next_ch, input_fasta_fh);

        if (next_ch == EOF)
        {
            break;
        }

        long beg = ftell(input_fasta_fh);
        fscanf(input_fasta_fh, "%*[^\n]\n");

        long end = ftell(input_fasta_fh);

        assert(beg < end);

        int retval = fseek(input_fasta_fh, beg, SEEK_SET);
        assert(retval == 0);

        if (next_ch == '>'){ // we have a header line
            ++chunk;
            if (chunk > 1)
            { 
                fprintf(output_fasta_fh, "\n");
            }
            //printing out a header line
            FileUtils::cat(line, MAX_LINE, end - beg, input_fasta_fh, output_fasta_fh);
        }
        else
        {
            //printing out a sequence chunk, without the newline
            FileUtils::cat(line, MAX_LINE, end - beg - 1, input_fasta_fh, output_fasta_fh);
            char nl = fgetc(input_fasta_fh); // get and throw away newline character
            assert(nl == '\n');
        }
    }

    fprintf(output_fasta_fh, "\n");

    fclose(input_fasta_fh);
    fclose(output_fasta_fh);

    delete line;
    
    return 0;
}

