#include "file_utils.h"

#include <cstring>

int main(int argc, char ** argv)
{
    char * infile = argv[1];
    size_t bytes_wanted = static_cast<size_t>(atof(argv[2]));
    size_t max_line_length = static_cast<size_t>(atof(argv[3]));

    FILE * input_fh = fopen(infile, "r");

    fseek(input_fh, 0L, SEEK_END);
    long end_pos = ftell(input_fh);
    fseek(input_fh, -1 * static_cast<long>(bytes_wanted), SEEK_END);
    
    char * buf = new char[max_line_length + 1];
    fgets(buf, max_line_length, input_fh);
    size_t bufsize = strlen(buf);
    if (buf[bufsize - 1] != '\n')
    {
        fprintf(stderr, "Error: Didn't find a newline within a max line length of %Zu bytes\n",
                max_line_length);
        exit(10);
    }
    delete buf;

    long cur_pos = ftell(input_fh);

    size_t catbuf_size = 100000000;
    char * catbuf = new char[catbuf_size + 1];
    
    FileUtils::cat(catbuf, catbuf_size, end_pos - cur_pos, input_fh, stdout);
    fclose(input_fh);
    delete catbuf;
    return 0;
}
