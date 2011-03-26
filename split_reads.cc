#include <cstring>
#include <vector>
#include <cassert>

#include "file_utils.h"
#include "dep/tools.h"
#include "cisortho/string_tools.h"
#include "zlib.h"


int print_usage(size_t mdef)
{
    fprintf(stderr, "\nUsage: split_reads [OPTIONS] split_list input.fq output_seg1.fq output_seg2.fq ...\n\n"
            "Options:\n"
            "-m  INT      number of bytes of memory to allocate [%Zu]\n"
            "'split_list' for example, '25,25,25' is the comma-separated list of read lengths to split input\n",
            mdef);
    return 1;
}

int main(int argc, char *argv[])
{
    size_t max_mem_def = 1024L * 1024L * 1024L; //1 GB ought to do it!
    size_t max_mem = max_mem_def;

    char c;
    while ((c = getopt(argc, argv, "m:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        default: return print_usage(max_mem_def); break;
        }
    }

    if (argc <= optind)
    {
        return print_usage(max_mem_def);
    }


    size_t const max_line = 4096;

    char const* split_list = argv[optind];
    std::vector<std::string> section_length_string = 
        split_token(split_list, strlen(split_list), ",");

    size_t num_sections = section_length_string.size();

    std::vector<size_t> section_length(num_sections);
    std::vector<size_t> section_offset(num_sections);
    section_offset[0] = 0;
    section_length[0] = atoi(section_length_string[0].c_str());
    for (size_t s = 1; s != num_sections; ++s)
    {
        section_length[s] = atoi(section_length_string[s].c_str());
        section_offset[s] = section_length[s-1] + section_offset[s-1];
    }

    int req_args = 2 + num_sections;
    int arg_count = optind + req_args;

    if (argc != arg_count)
    {
        return print_usage(max_mem_def);
    }

    char const* input_fq_file = argv[optind + 1];
    char *const* output_fq_files = argv + optind + 2;

    FILE * input_fq_fh = open_or_die(input_fq_file, "r", "Fastq input file");

    FILE ** output_fq_fh = new FILE *[num_sections];

    gzFile fastq_zfh = gzopen(input_fq_file, "r");
    std::vector<size_t> chunk_lengths = 
        FileUtils::ChunkLengths(fastq_zfh, max_mem, max_line);
    gzclose(fastq_zfh);
    
    char * chunk_buffer = new char[max_mem + max_line];
    char * out_chunk_buffer = new char[max_mem + max_line];

    size_t max_section_size = max_mem / num_sections;

    char ** out_ptr_start = new char *[num_sections];
    char ** out_ptr_end = new char *[num_sections];
    char ** out_ptr = new char *[num_sections];

    for (size_t s = 0; s != num_sections; ++s)
    {
        out_ptr_start[s] = out_chunk_buffer + (s * max_section_size);
        out_ptr_end[s] = out_ptr_start[s] + max_section_size - max_line;
    }

    size_t num_chunks = chunk_lengths.size();
    size_t readlength;

    for (size_t s = 0; s != num_sections; ++s)
    {
        output_fq_fh[s] = open_or_die(output_fq_files[s], "w", "Output split fastq file");
    }

    std::vector<char *> fastq_lines;

    char * left_over_part = new char[max_line];
    left_over_part[0] = '\0';

    for (size_t chunk = 0; chunk != num_chunks; ++chunk)
    {
        strcpy(chunk_buffer, left_over_part);
        size_t left_over_bytes = strlen(left_over_part);
        
        size_t nbytes_read = fread(chunk_buffer + left_over_bytes, 1, 
                                   chunk_lengths[chunk], input_fq_fh);

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
        
        char * read;
        char * id;
        char idstring[256];
        char * qualstring;
        float f_readlength;

        for (size_t s = 0; s != num_sections; ++s)
        {
            out_ptr[s] = out_ptr_start[s];
        }

        std::vector<char *>::iterator fit;
        for (fit = fastq_lines.begin(); fit != chunk_end; ++fit)
        {
            id = *fit++; ++id;
            read = *fit++; ++fit;
            readlength = strlen(read);
            f_readlength = static_cast<float>(readlength);
            qualstring = *fit;

            //write each section of this read to memory buffer
            for (size_t s = 0; s != num_sections; ++s)
            {
                char *& ptr = out_ptr[s];
                sprintf(idstring, "@%s|%Zu:%Zu\n", id, section_offset[s], s);
                strcpy(ptr, idstring);
                ptr += strlen(idstring);
                strncpy(ptr, read + section_offset[s], section_length[s]);
                ptr += section_length[s];
                strcpy(ptr, "\n+\n");
                ptr += 3;
                strncpy(ptr, qualstring + section_offset[s], section_length[s]);
                ptr += section_length[s];
                strcpy(ptr, "\n");
                ++ptr;

                //flush each memory buffer as necessary
                if (out_ptr_end[s] < ptr)
                {
                    size_t num_chars = std::distance(out_ptr_start[s], ptr);
                    fwrite(out_ptr_start[s], 1, num_chars, output_fq_fh[s]);
                    fflush(output_fq_fh[s]);
                    out_ptr[s] = out_ptr_start[s];
                }
            }
        }
        for (size_t s = 0; s != num_sections; ++s)
        {
            //flush each memory buffer as necessary
            fwrite(out_ptr_start[s], 1, std::distance(out_ptr_start[s], out_ptr[s]),
                   output_fq_fh[s]);
            
            out_ptr[s] = out_ptr_start[s];
        }
        
    }

    for (size_t s = 0; s != num_sections; ++s)
    {
        fclose(output_fq_fh[s]);
    }

    delete chunk_buffer;
    delete out_chunk_buffer;
    delete left_over_part;
    delete out_ptr;
    delete out_ptr_start;
    delete out_ptr_end;
    delete output_fq_fh;

	return 0;
}
