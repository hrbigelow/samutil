#include "sam_helper.h"

int usage()
{
    fprintf(stderr, "Usage: samutil_check input.sam\n");
    return 1;
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        return usage();
    }

    char const* sam_file = argv[1];
    FILE * sam_fh = fopen(sam_file, "r");

    bool source_is_ones_based = true;
    SamLine * samline;
    while ((samline = new SamLine(sam_fh, source_is_ones_based)))
    {
        if (samline->parse_flag != HEADER)
        {
            break;
        }
        delete samline;
    }

    while (samline->parse_flag == DATA_LINE)
    {
        delete samline;
        samline = new SamLine(sam_fh, source_is_ones_based);
    }
    if (samline->parse_flag == END_OF_FILE)
    {
        fprintf(stdout, "Passed\n");
    }
    else
    {
        fprintf(stdout, "Failed\n");
    }
    return 0;
}
