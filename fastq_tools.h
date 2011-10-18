//determine whether fastq file is a 33 or 64 offset, store in qual_offset
//if scale is invalid, set qual_offset to -1 and return false
#include <cstdio>

bool fastq_file_offset(FILE * fastq_fh, int * qual_offset, char ** file_types);

// samples first, mid, and last qual chars of every line of fastq_fh, saving them
// in minc and maxc
void fastq_extreme_chars(FILE * fastq_fh, char * minc, char * maxc);
