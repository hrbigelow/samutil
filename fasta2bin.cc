#include <cstdio>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "dep/bindepth.h"

/*
  fasta bin format is:

  file: contig-dictionary seq-data [, seq-data ...]

  contig-dictionary: (see ~/cc/dep/bindepth.cc)

  seq-data: character array holding all sequence (no newlines)

 */

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        fprintf(stderr, "Converts a fasta file to binary format for easy processing\n"
                "Usage: fasta2bin contig_dict.rdb in.fasta out.fbin\n");
        return 1;
    }
    const char *contig_dict_file = argv[1], *in_file = argv[2], *out_file = argv[3];
    
    FILE *contig_dict_fh, *in_fh, *out_fh;
    if (! (contig_dict_fh = fopen(contig_dict_file, "r")))
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", contig_dict_file);
        return 1;
    }
    if (! (in_fh = fopen(in_file, "r")))
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", in_file);
        return 1;
    }
    if (! (out_fh = fopen(out_file, "w")))
    {
        fprintf(stderr, "Error: couldn't open output file %s\n", out_file);
        return 1;
    }

    size_t ncontigs;
    contig_dict_t *contigs = parse_contig_dict(contig_dict_fh, &ncontigs), *cp = contigs;
    fclose(contig_dict_fh);

    write_contig_dict(contigs, ncontigs, out_fh);

    int ch;
    char contig[20];
    // prescan whole file
    while (1)
    {
        ch = fgetc(in_fh); assert((char)ch == '>');
        fgets(contig, 20, in_fh);
        contig[strlen(contig) - 1] = '\0'; // replace \n with \0
        if (strcmp(contig, cp->name))
        {
            fprintf(stderr, "Error: next contig in fasta file is %s, but contig dictionary has %s.  Breaking.\n",
                    cp->name, contig);
            break;
        }
        // now write out everything up until EOF or '>'
        while ((char)(ch = fgetc(in_fh)) != '>' && ch != EOF)
        {
            if (ch != '\n') fputc(ch, out_fh);
        }
        if (ch == EOF) break;

        ungetc(ch, in_fh);
        cp++;
    }

    free((void *)contigs);
    fclose(in_fh);
    fclose(out_fh);
    return 0;
}
