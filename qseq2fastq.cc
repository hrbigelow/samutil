//Convert Illumina GAPipeline 1.3.2 qseq format to FastQ format, as a filter
/*
  Compile with:
  g++ qseq2fastq.cpp -o qseq2fastq
  Use as:
  /path/to/qseq2fastq < s_1_qseq.txt > s_1.fastq
*/
#include <algorithm>

int main(int argc, char ** argv){

    if (argc == 0)
    {
        printf("usage: qseq2fastq < input_qseq.txt > output.fastq\n");
        return 0;
    }

    int lane, tile, x, y, read;
    char machine[256];
    char sequence[1024];
    char qseq_string[1024];
    char fastq_string[1024];

    char qseq2fastq[256];
//     for (int i=0; i < 256; ++i){
//         qseq2fastq[i] = char(std::max(0,i-31));
//     }

    char *q, *f;
    int nfields_read;
    while (! feof(stdin)){
        nfields_read = 
        fscanf(stdin, "%s %*s %i %i %i %i %*i %i %s %s %*i\n", machine, &lane, &tile, &x, &y, &read, sequence, qseq_string);
        if (nfields_read != 8)
        {
            if (nfields_read > 0)
            {
                fprintf(stderr, "qseq2fastq: encountered badly formatted input with %i fields.\n", nfields_read);
                exit(1);
            }
            break;
        }
        //for (q = qseq_string, f = fastq_string; *q != 0; ++q, ++f){ *f = qseq2fastq[int(*q)]; }
        //*f = 0;
        fprintf(stdout, "@%s:%i:%i:%i:%i/%i\n%s\n+\n%s\n", machine, lane, tile, x, y, read, sequence, qseq_string);
    }
    return 0;
}
