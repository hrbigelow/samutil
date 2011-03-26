#include "sim.h"

#include <cstdio>
#include <cstring>

int sim_usage()
{
    fprintf(stderr,
            "\nAuthor: Henry Bigelow (hbigelow@amgen.com)\n\n"
            "Usage:\n\n"
            "sim reads        Simulate paired reads\n"
            "sim expression   Simulate transcript expression levels\n\n"
            );
    return 1;
}


int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        return sim_usage();
    }
    else if (strcmp(argv[1], "reads") == 0)
    {
        return main_sim_reads(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "expression") == 0)
    {
        return main_sim_expression(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
}
