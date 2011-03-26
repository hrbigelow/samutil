#include <cstring>
#include <vector>
#include <cassert>

#include "file_utils.h"
#include "fastq_tools.h"
#include "dep/tools.h"
#include "zlib.h"


int const base_to_index[] =
    {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
    };
char const* bases_upper = "ACGTN";



void format_qual_string(char const* qual, char * out_qual, int qual_value_adjust)
{
    size_t l = strlen(qual);
    if (qual_value_adjust == 0)
    {
        strcpy(out_qual, qual);
    }
    else
    {
        for (size_t i = 0; i != l; ++i)
        {
            out_qual[i] = qual[i] + qual_value_adjust;
        }
    }
}


size_t filter_garbage_reads(char const* fastq_file, 
                            char const* filtered_fastq_file,
                            size_t max_mem, size_t max_line,
                            float max_frac_base[5],
                            int qual_value_adjust)
{	
    gzFile fastq_zfh = gzopen(fastq_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(fastq_zfh, max_mem, max_line);
    gzclose(fastq_zfh);
    
    char * chunk_buffer = new char[max_mem + max_line];
    char * out_chunk_buffer = new char[max_mem + max_line];
    char * out_ptr_end = out_chunk_buffer + max_mem;
    char * out_ptr;

    size_t num_chunks = chunk_lengths.size();
    char formatted_qual[max_line];
    size_t next_id = 0;
    size_t num_reads = 0;
    size_t num_reads_chucked = 0;
    size_t readlength;


    FILE * fastq_fh = open_or_die(fastq_file, "r", "Input fastq file");
    FILE * filtered_fastq_fh = open_or_die(filtered_fastq_file, "w", "Output filtered fastq file");

    std::vector<char *> fastq_lines;

    char * left_over_part = new char[max_line];
    left_over_part[0] = '\0';

    for (size_t chunk = 0; chunk != num_chunks; ++chunk)
    {
        out_ptr = out_chunk_buffer;
        strcpy(chunk_buffer, left_over_part);
        size_t left_over_bytes = strlen(left_over_part);

        size_t nbytes_read = fread(chunk_buffer + left_over_bytes, 1, 
                                   chunk_lengths[chunk], fastq_fh);

        assert(nbytes_read == chunk_lengths[chunk]);
        chunk_buffer[nbytes_read + left_over_bytes] = '\0';

        fastq_lines = FileUtils::find_line_starts(chunk_buffer);

        char * line_end = chunk_buffer;
        while ((line_end = strchr(line_end, '\n')) != NULL)
        {
            *line_end++ = '\0';
        }

        size_t num_left_over = fastq_lines.size() % 4;

        std::vector<char *>::iterator chunk_end = fastq_lines.end() - num_left_over;
        std::vector<char *>::iterator chunk_iter;

        for (chunk_iter = chunk_end; chunk_iter != fastq_lines.end(); ++chunk_iter)
        {
            line_end = *chunk_iter + strlen(*chunk_iter);
            *line_end = '\n';
        }
        char const* terminal = (chunk_end == fastq_lines.end()) ? "" : *chunk_end;
        strcpy(left_over_part, terminal);

        float counts[5];
        char * read;
        char * read_ptr;
        char * id;
        char * qualstring;
        float f_readlength;

        std::vector<char *>::iterator fit;
        for (fit = fastq_lines.begin(); fit != chunk_end; ++fit)
        {
            id = *fit++; ++id;
            read = *fit++; ++fit;
            readlength = strlen(read);
            f_readlength = static_cast<float>(readlength);
            qualstring = *fit;

            ++num_reads;
            ++next_id;
            
            std::fill(counts, counts + 5, 0.0);
            for (read_ptr = read; *read_ptr != '\0'; ++read_ptr)
            {
                counts[base_to_index[static_cast<int>(*read_ptr)]]++;
            }
            if (counts[0] / readlength > max_frac_base[0] ||
                counts[1] / readlength > max_frac_base[1] ||
                counts[2] / readlength > max_frac_base[2] ||
                counts[3] / readlength > max_frac_base[3] ||
                counts[4] / readlength >= max_frac_base[4])
            {
                ++num_reads_chucked;
                continue;
            }

            for (read_ptr = read; *read_ptr != '\0'; ++read_ptr)
            {
                *read_ptr = bases_upper[base_to_index[static_cast<int>(*read_ptr)]];
            }
            
            format_qual_string(qualstring, formatted_qual, qual_value_adjust);

            if (out_ptr_end < out_ptr)
            {
                //flush output, reset pointer
                fwrite(out_chunk_buffer, 1, 
                       std::distance(out_chunk_buffer, out_ptr),
                       filtered_fastq_fh);
                out_ptr = out_chunk_buffer;
            }
            int num_printed = sprintf(out_ptr, "@%Zu\n%s\n+%s\n%s\n", 
                                      next_id, read, id, formatted_qual);
            out_ptr += num_printed;
        }

        //flush remaining output
        fwrite(out_chunk_buffer, 1, 
               std::distance(out_chunk_buffer, out_ptr),
               filtered_fastq_fh);
        out_ptr = out_chunk_buffer;
    }

    delete chunk_buffer;
    delete left_over_part;

    fclose(fastq_fh);
    fclose(filtered_fastq_fh);

    return num_reads_chucked;
}


int print_usage(char const* fdef, int qdef, size_t mdef)
{
    fprintf(stderr, 
            "Usage: filter_reads_for_tophat [OPTIONS] input.fq filtered_output.fq\n"
            "Options:\n\n"
            "-f  STRING   base composition filter string.  e.g. 0.9,0.9,0.9,0.9,0.1\n"
            "             providing maximum fraction of each base (A,C,G,T,N)\n"
            "             for a read to be kept [%s]\n"
            "-q  INT      desired output quality score offset [%d]\n"
            "-m  INT      number of bytes of memory to allocate [%Zu]\n",
            fdef, qdef, mdef);
    return 1;
}

int main(int argc, char *argv[])
{
    char const* max_frac_base_string_def = "0.9,0.9,0.9,0.9,0.1";
    int output_qual_offset_def = 33;
    char const* max_frac_base_string = max_frac_base_string_def;
    int output_qual_offset = output_qual_offset_def;

    size_t max_mem_def = 1024L * 1024L * 1024L; //1 GB ought to do it!
    size_t max_mem = max_mem_def;

    char c;
    while ((c = getopt(argc, argv, "f:q:m:")) >= 0)
    {
        switch(c)
        {
        case 'f': max_frac_base_string = optarg; break;
        case 'q': output_qual_offset = static_cast<int>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        default: return print_usage(max_frac_base_string_def, 
                                    output_qual_offset_def,
                                    max_mem_def); break;
        }
    }

    int req_args = 2;
    int arg_count = optind + req_args;

    if (argc != arg_count)
    {
        return print_usage(max_frac_base_string_def, 
                           output_qual_offset_def, 
                           max_mem_def);
    }

    float max_frac_base[5];
    int num_parsed = sscanf(max_frac_base_string, "%f,%f,%f,%f,%f", 
                            max_frac_base,
                            max_frac_base + 1,
                            max_frac_base + 2,
                            max_frac_base + 3,
                            max_frac_base + 4);

    if (num_parsed != 5)
    {
        fprintf(stderr, "Error: option -f (%s) has bad format.\n"
                "Should have five comma-separated numbers\n",
                max_frac_base_string);
        exit(1);
    }

    char const* fastq_file = argv[optind];
    char const* filtered_fastq_file = argv[optind + 1];

    FILE * fastq_fh = open_or_die(fastq_file, "r", "Fastq input file");

    size_t max_line = 4096L * 4;

    int input_qual_offset;
    bool fastq_is_valid = fastq_file_offset(fastq_fh, &input_qual_offset);
    if (! fastq_is_valid)
    {
        fprintf(stderr, "Error: couldn't determine fastq file %s quality scale\n", fastq_file);
        exit(1);
    }

    fclose(fastq_fh);
    int qual_value_adjust = output_qual_offset - input_qual_offset;

    size_t num_reads_chucked =
        filter_garbage_reads(fastq_file, filtered_fastq_file, 
                             max_mem, max_line,
                             max_frac_base, qual_value_adjust);
    
    fprintf(stderr, "Filtered out %Zu reads\n", num_reads_chucked);

	return 0;
}
