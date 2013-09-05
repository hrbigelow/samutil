#include <cstdio>
#include <cstring>
#include <vector>

#include "file_utils.h"
#include "dep/tools.h"

int deal_fastq_usage(size_t mdef)
{
    fprintf(stderr, 
            "\nUsage:\n\n"
            "deal_fastq [OPTIONS] in.fq out1.fq [out2.fq [...]]\n\n"
            "De-intercalates ('deals') a fastq file. The user must know\n"
            "How many matching reads there are\n\n"
            "Options:\n\n"
            "-m     INT     number bytes of memory to use [%zu]\n\n",
            mdef);
    return 1;
}

int main(int argc, char ** argv)
{
    size_t mdef = 1024l * 1024l * 256l; // 256 MB memory
    size_t max_mem = mdef;

    char c;
    while ((c = getopt(argc, argv, "m:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        default: return deal_fastq_usage(mdef); break; 
        }
    }

    if (argc - optind < 2)
    {
        return deal_fastq_usage(mdef);
    }

    FILE * in_fh = open_or_die(argv[optind], "r", "Input fastq file");
    size_t n_outfiles = argc - optind - 1;

    size_t in_mem = max_mem * 4 / 10; // 40% for in file buffer
    size_t out_mem = (max_mem * 6) / (n_outfiles * 10); // (60 % divided among out file buffers
    // this uneven split is to guard against possible differently lengthed paired reads.

    char * chunk_buffer_in = new char[in_mem + 1];
    char ** chunk_buffer_out = new char*[n_outfiles];
    char ** write_pointer = new char*[n_outfiles];
    
    FILE ** out_fhs = new FILE*[n_outfiles];
    for (size_t o = 0; o != n_outfiles; ++o)
    {
        out_fhs[o] = open_or_die(argv[optind + o + 1], "w", "Output fastq file");
        chunk_buffer_out[o] = new char[out_mem + 1];
    }

    char * last_fragment;
    char * read_pointer = chunk_buffer_in;
    size_t nbytes_read, nbytes_unused = 0;

    std::vector<char *> lines;
    std::vector<char *>::const_iterator lit, lit2;

    while (! feof(in_fh))
    {
        nbytes_read = fread(read_pointer, 1, in_mem - nbytes_unused, in_fh);
        read_pointer[nbytes_read] = '\0';
        lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        size_t n_remaining = lines.size() % (4 * n_outfiles);

        for (size_t o = 0; o != n_outfiles; ++o)
        {
            write_pointer[o] = chunk_buffer_out[o];
            size_t const incr = 4 * n_outfiles;
            size_t offset = o * 4;
            for (lit = lines.begin(); lit != lines.end() - n_remaining; lit += incr)
            {
                lit2 = lit + offset;
                assert((*lit2)[0] == '@');
                write_pointer[o] += sprintf(write_pointer[o], "%s\n%s\n+\n%s\n",
                                            *(lit2), *(lit2 + 1), *(lit2 + 3));
            }

            fputs(chunk_buffer_out[o], out_fhs[o]);
        }
        
        // recycle remaining lines.  the 'end' pointer also represents
        // the end of the buffer
        char * unused_start = lit == lines.end() ? last_fragment : *(lit);
        
        for ( ; lit != lines.end(); ++lit)
        {
            assert(*((*lit) + strlen(*lit)) == '\0');
            *((*lit) + strlen(*lit)) = '\n';
        }
        nbytes_unused = strlen(unused_start);

        memmove(chunk_buffer_in, unused_start, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
    }
        

    for (size_t c = 0; c != n_outfiles; ++c)
    {
        fclose(out_fhs[c]);
        delete chunk_buffer_out[c];
    }
    delete [] out_fhs;
    delete [] chunk_buffer_in;
    delete [] chunk_buffer_out;
}
