#include "file_utils.h"

int usage()
{
    fprintf(stderr, 
            "\nUsage:\n"
            "file_bytes get <input_file> <chunk_size> <chunk_num>\n\n"
            "returns the text of the chunk chosen (chunk numbered from zero)\n\n"
            "file_bytes count <input_file> <chunk_size>\n"
            "returns the number of chunks\n"
            );
    return 1;
}

int main(int argc, char ** argv)
{
    if (argc < 4)
    {
        return usage();
    }

    char * mode = argv[1];
    char * input_file = argv[2];
    size_t chunk_bytes = static_cast<size_t>(atof(argv[3]));


    if (strcmp(mode, "count") == 0)
    {
        if (argc != 4)
        {
            return usage();
        }
        FILE * input_fh = fopen(input_file, "r");
        fprintf(stdout, "%Zu\n", count_chunks(input_fh, chunk_bytes));
        fclose(input_fh);
    }
    else if (strcmp(mode, "get") == 0)
    {
        if (argc != 5)
        {
            return usage();
        }
        FILE * input_fh = fopen(input_file, "r");
        size_t chunk_num = static_cast<size_t>(atof(argv[4]));
        print_chunk_by_line(input_fh, chunk_num, chunk_bytes, stdout);
        fclose(input_fh);
    }
    else
    {
        return usage();
    }
    return 0;
}

