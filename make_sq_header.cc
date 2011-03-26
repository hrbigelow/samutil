#include <cstdio>
#include <cstdlib>
#include <fstream>

//produce a '@SQ\tSN:...\tLN:# header from a fasta file
int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "make_sq_header contigs.fa header.sam\n");
        return 1;
    }

    char * fasta = argv[1];
    char * header = argv[2];
    FILE * fasta_fh = fopen(fasta, "r");
    if (fasta_fh == NULL)
    {
        fprintf(stderr, "Couldn't open fasta file %s\n", fasta);
        exit(1);
    }
    FILE * header_fh = fopen(header, "w");

    char contig[1000];
    size_t contig_start;
    size_t contig_end;
    size_t contig_cur;
    size_t contig_length;
    size_t num_newlines;

    while (! feof(fasta_fh))
    {
        int nfields = fscanf(fasta_fh, ">%s\n", contig);
        if (nfields != 1)
        {
            fprintf(stderr, "Error: bad fasta format\n");
            exit(1);
        }

        contig_start = ftell(fasta_fh);
        fscanf(fasta_fh, "%*[^>]");
        contig_end = ftell(fasta_fh); //includes newline
        fseek(fasta_fh, contig_start, std::ios::beg);
        contig_cur = contig_start;
        num_newlines = 0;
        while (contig_cur != contig_end)
        {
            fscanf(fasta_fh, "%*[^\n]\n");
            contig_cur = ftell(fasta_fh);
            ++num_newlines;
        }
        contig_length = contig_end - contig_start - num_newlines;
        fprintf(header_fh, "@SQ\tSN:%s\tLN:%Zu\n", contig, contig_length);
    }
    
    fclose(fasta_fh);
    fclose(header_fh);
}
