#include <cstdlib>
#include <cstdio>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

int main(int argc, char **argv){

  int truncate_depth = atoi(argv[1]);
  int depth;
  int64_t * counts = new int64_t[truncate_depth+1];
  for (int i = 0; i <= truncate_depth; ++i){ counts[i] = 0; }
  while (! feof(stdin)){
    scanf("%i\n", &depth);
    counts[depth > truncate_depth ? truncate_depth : depth]++;
  }

  for (int i = 0; i <= truncate_depth; ++i){ 
    printf("%i\t%"PRId64"\n", i, counts[i]);
  }

  delete counts;
}

