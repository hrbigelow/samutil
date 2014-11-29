#include "sam_truncate.h"
#include "file_utils.h"
#include "sam_line.h"
#include "dep/tools.h"

#include <vector>
#include <cstdio>
#include <unistd.h>
#include <omp.h>
#include <parallel/algorithm>

int sam_truncate_usage(size_t mdef)
{
    fprintf(stderr,
            "\nUsage:\n\n"
            "samutil truncate [OPTIONS] in.sam out.sam\n\n"
            "Options:\n\n"
            "-t     INT     number of threads to use [1]\n"
            "-m     INT     number bytes of memory to use [%Zu]\n"
            "-C     STRING  work in the directory named here [.]\n"
            "-a     FLAG    if present, append to output file. (otherwise create fresh) [false]\n"
            "\n"
            "Notes:\n\n"
            "Abbreviates QNAME, SEQ and QUAL. Converts CIGAR '*' to 'nS'\n\n",
            mdef
            );
    return 1;
}


// input: full_samline raw string
// output: in-place modified string containing truncated string
struct truncate_sam_unary
{
    truncate_sam_unary() { }
    char * operator()(char * full_samline)
    {
        SamLine * sam = new SamLine(full_samline);
        if (sam->parse_flag == PARSE_ERROR)
        {
            fprintf(stderr, "Encountered error parsing input\n");
            exit(1);
        }
        
        sam->sprint(full_samline);
        delete sam;
        return full_samline;
    }
};







int main_sam_truncate(int argc, char ** argv)
{

    size_t max_mem_def = 1024l * 1024l * 1024l * 4l; // 4 GB memory
    size_t max_mem = max_mem_def;

    size_t num_threads = 1;
    char const* working_dir = ".";

    bool append_output = false;

    char c;
    while ((c = getopt(argc, argv, "m:t:C:a")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'C': working_dir = optarg; break;
        case 'a': append_output = true; break;
        default: return sam_truncate_usage(max_mem_def); break;
        }
    }

    int arg_count = optind + 2;
    if (argc != arg_count)
    {
        return sam_truncate_usage(max_mem_def);
    }

    char const* input_sam_file = argv[optind];
    char const* output_sam_file = argv[optind + 1];

    int chdir_success = chdir(working_dir);
    if (chdir_success != 0)
    {
        fprintf(stderr, "Error: couldn't change directory to %s\n", working_dir);
        exit(1);
    }

    FILE * input_sam_fh = open_or_die(input_sam_file, "r", "Input SAM file to be truncated");

    char const* out_mode = append_output ? "a" : "w";

    FILE * output_sam_fh = open_or_die(output_sam_file, out_mode, "Truncated output SAM file");

    omp_set_dynamic(false);
    omp_set_num_threads(num_threads);


    char * header_buf = ReadAllocSAMHeader(input_sam_fh);
    size_t header_length = strlen(header_buf);
    fwrite(header_buf, 1, header_length, output_sam_fh);
    delete [] header_buf;
    
    size_t chunk_size = max_mem / 2;
    char * chunk_buffer_in = new char[chunk_size + 1];
    char * chunk_buffer_out = new char[chunk_size + 1];

    size_t nbytes_unused = 0;
    size_t nbytes_read;
    std::vector<char *> sam_lines;
    std::vector<char *>::iterator sit, end;

    char * write_pointer;
    char * read_pointer = chunk_buffer_in;
    char * last_fragment;

    while (! feof(input_sam_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, input_sam_fh);
        read_pointer[nbytes_read] = '\0';
        sam_lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        if (! SamLine::initialized && ! sam_lines.empty())
        {
            char const* dummy_read_layout = "y"; // merely to tell SamLine to parse traditional SAM.
            SamLine::SetGlobalFlags(QNAMEFormat(sam_lines[0]), dummy_read_layout, "", 0, false);
        }

        // do not use the last line. it is incomplete.
        truncate_sam_unary truncate_aux;
        __gnu_parallel::transform(sam_lines.begin(), sam_lines.end(), 
                                  sam_lines.begin(), truncate_aux);

        write_pointer = chunk_buffer_out;
        for (sit = sam_lines.begin(); sit != sam_lines.end(); ++sit)
        {
            strcpy(write_pointer, (*sit));
            write_pointer += strlen(*sit);
        }
        fwrite(chunk_buffer_out, 1, write_pointer - chunk_buffer_out, output_sam_fh);

        // now restore unused bytes
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

    }
    delete chunk_buffer_in;
    delete chunk_buffer_out;

    fclose(input_sam_fh);
    fclose(output_sam_fh);
    
    return 0;
}
