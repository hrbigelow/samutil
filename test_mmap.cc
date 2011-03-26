#include <sys/mman.h>
#include <cstdlib>
#include <sys/timeb.h>

int main(int argc, char ** argv)
{
    
    timeb millitime;
    ftime(& millitime);

    size_t N = static_cast<size_t>(atof(argv[0]));
    int seed = 98398;

    fprintf(stdout, "time 1: %i\n", millitime.millitm);
    int * rands = new int[N];
    srand(seed);
    for (size_t i = 0; i != N; ++i)
    {
        rands[i] = rand();
    }
    std::sort(rands, rands + N);
    delete rands;

    ftime(& millitime);
    fprintf(stdout, "time 2: %i\n", millitime.millitm);

    char * tmp_file_template = "rands.XXXXXX";

    int file_des = mkstemp(tmp_file_template);
    FILE * fh  = fdopen(file_des, "r+");

    void * vbuf = mmap(NULL, chunk_size, PROT_READ | PROT_WRITE, 
                       MAP_PRIVATE, file_des, 0);

    int * mrands = vbuf;
    srand(seed);
    for (size_t i = 0; i != N; ++i)
    {
        fwrite(mrands[i] = rand();
    }
    std::sort(mrands, mrands + N);

    ftime(& millitime);
    fprintf(stdout, "time 3: %i\n", millitime.millitm);

    return 0;
}
