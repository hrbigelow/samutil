#include "sam_file.h"

#include <vector>

// reads and allocates a buffer loaded with the sam header
char * ReadAllocSAMHeader(FILE * sam_fh)
{
    std::vector<char> header;
    header.reserve(1000);

    char p,n;
    p = '\0';

    while (! feof(sam_fh))
    {
        n = fgetc(sam_fh);
        if ((p == '\0' || p == '\n') && n != '@')
        {
            ungetc(n, sam_fh);
            clearerr(sam_fh);
            break;
        }
        p = n;
        header.push_back(p);
    }
    char * header_string = new char[header.size() + 1];
    std::copy(header.begin(), header.end(), header_string);
    header_string[header.size()] = '\0';

    return header_string;
}


//sets sam_fh to position of first data line
void SetToFirstDataLine(FILE ** sam_fh)
{
    rewind(*sam_fh);

    //set just after first newline that doesn't begin with '@'
    char c;
    size_t dummy;
    char *line = NULL;

    while ((c = fgetc(*sam_fh)) == '@')
    {
        ungetc(c, *sam_fh);
        getline(&line, &dummy, *sam_fh);
    }
    ungetc(c, *sam_fh);

    return;
}


void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh)
{
    SetToFirstDataLine(input_sam_fh);

    size_t header_length = ftell(*input_sam_fh);
    rewind(*input_sam_fh);

    char * header_buf = new char[header_length];
    fread(header_buf, 1, header_length, *input_sam_fh);
    fwrite(header_buf, 1, header_length, output_fh);
    delete [] header_buf;
    fflush(output_fh);
}
