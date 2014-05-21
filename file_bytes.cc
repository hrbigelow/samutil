#include "file_utils.h"

#include <cstring>

int usage()
{
    fprintf(stderr, 
            "\nUsage:\n"
            "file_bytes get <input_file> <chunk_size> <max_line_length> <chunk_number>\n\n"
            "returns the text of the chunk chosen (chunk numbered from zero)\n\n"
            "file_bytes count <input_file> <chunk_size> <max_line_length>\n"
            "returns the number of chunks\n"
            );
    return 1;
}

int main(int argc, char ** argv)
{
    if (argc < 5)
    {
        return usage();
    }

    char * mode = argv[1];
    char * input_file = argv[2];
    size_t chunk_bytes = static_cast<size_t>(atof(argv[3]));
    size_t max_line_length = static_cast<size_t>(atof(argv[4]));

    if (strcmp(mode, "count") == 0)
    {
        if (argc != 5)
        {
            return usage();
        }
        FILE * input_fh = fopen(input_file, "r");
        std::vector<size_t> chunk_lengths = FileUtils::chunk_lengths(input_fh, chunk_bytes, max_line_length);
        fprintf(stdout, "%Zu\n", chunk_lengths.size());
        fclose(input_fh);
    }
    else if (strcmp(mode, "get") == 0)
    {
        if (argc != 6)
        {
            return usage();
        }
        FILE * input_fh = fopen(input_file, "r");
        size_t chunk_num = static_cast<size_t>(atof(argv[5]));
        std::vector<size_t> chunk_lengths = FileUtils::chunk_lengths(input_fh, chunk_bytes, max_line_length);
        long offset = 0;
        for (size_t c = 0; c != chunk_num; ++c)
        {
            offset += chunk_lengths[c];
        }
        fseek(input_fh, offset, SEEK_SET);
        size_t buf_size = 1000000;
        char * buf = new char[buf_size + 1];
        FileUtils::cat(buf, buf_size, chunk_lengths[chunk_num], input_fh, stdout);
        fclose(input_fh);
    }
    else
    {
        return usage();
    }
    return 0;
}

