#include <cstdarg>
#include <cstdio>

void Variadic(char const* format,
              char const* line,
              ...)
{
    va_list args;
    va_start(args, line);
    int count = vsscanf(line, format, args);
    va_end(args);
    
}


int main(int argc, char ** argv)
{
    char line[] = "field1\tfield2\n";
    char field1[100];
    char field2[100];
    void * vfield1 = field1;
    void * vfield2 = field2;

    Variadic("%i\t%i\n", line, vfield1, vfield2);
    
    return 0;
}
