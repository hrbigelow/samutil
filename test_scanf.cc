#include <cstdio>

int main()
{
    char * pat = "a:b:23:d\n";
    char a, b;
    int i;
    sscanf(pat, "%c:%c:%i", &a, &b, &i);
    printf("%c\t%c\t%i\n", a, b, i);
    return 0;
}
