//output character counts of maq pileup string containing "acgtACGT,."
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

const int MAX_PILEUP_DEPTH=10000000;

int sum_of_reads(char const* accepted, int const* char_index, 
                 int const* tally){
  int sum = 0;
  for (char const* c = accepted; *c != 0; ++c) {
    int ci = static_cast<int>(*c);
    sum += tally[char_index[ci]];
  }
  return sum;
}


int main (int argc, char ** argv){

  char const* chars_to_index = argv[1]; // i.e. 'ACGTacgt.,'
  int min_non_wt_reads = atoi(argv[2]);

  int char_index[256];
  int num_chars = std::strlen(chars_to_index);

  for (int c=0; c < 256; ++c) char_index[c] = 255;
  for (int c=0; c < num_chars; ++c)
    char_index[static_cast<int>(chars_to_index[c])] = c;

  char * pileup = new char[MAX_PILEUP_DEPTH];
  char * ptr;
  int tally[256];

  char id[1000];
  int64_t position;
  char snp;
  int depth;
  
  while (! feof(stdin)){

    scanf("%s\t%"PRId64"\t%c\t%i\t%s\n", id, &position, &snp, &depth, pileup);
    ptr = pileup;
    //printf("%s\n", pileup);

    for (int t=0; t < num_chars; ++t) tally[t] = 0;

    while (*ptr != 0) {
      tally[char_index[static_cast<int>(*ptr++)]]++;
      //printf("%c\n", *ptr);
    }

//     printf("non-wt reads: %i, wt-reads: %i\n", 
//             sum_of_reads("ACGTacgt", char_index, tally),
//             sum_of_reads(",.", char_index, tally));

    if (sum_of_reads("ACGTacgt", char_index, tally) < min_non_wt_reads) continue;

    printf("%s\t%"PRId64"\t%c\t%i\t%s", id, position, snp, depth, pileup);
    
    for (int t=0; t < num_chars; ++t)
      printf("\t%i", tally[t]);

    printf("\n");
  }

  delete pileup;

  return 0;

}
