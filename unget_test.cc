#include <cstdio>


int main(int argc, char ** argv)
{
    FILE * fh = fopen(argv[1], "r");

    char * buf = new char[1000];
    while (! feof(fh))
    {
        int nbytes_read = fread(buf, 1, 1000, fh);
        
        printf("feof flag is %s set after %i bytes read\n", (feof(fh) ? "indeed" : "not"), nbytes_read);

        // char test_char = getc(fh);
        // ungetc(test_char, fh);
        // bool is_last_chunk = feof(fh);
        
        // printf("after char test, feof flag is %s set\n", (is_last_chunk ? "indeed" : "not"));

    }

    delete buf;
    return 0;

}
