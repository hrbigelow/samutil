#include <stdio.h>
#include <unistd.h>


// if file is not NULL and not '/dev/null', attempts to
// open the file.  In this case, it is an error if it cannot open the file.
// otherwise, returns NULL
FILE * open_if_present(const char *file, const char *mode)
{
    if (file == NULL
        || strcmp(file, "/dev/null") == 0
        || strcmp(file, "") == 0)
    {
        return NULL;
    }
    else
    {
        FILE * fh = fopen(file, mode);
        if (fh == NULL)
        {
            fprintf(stderr, "Error: open_if_present: file %s not blank or '/dev/null' but"
                    " still couldn't open it\n", file);
            exit(1);
        }
        else
        {
            return fh;
        }
    }
}


int close_if_present(FILE *fh)
{
    if (fh != NULL) return fclose(fh);
    else return 0;
}
