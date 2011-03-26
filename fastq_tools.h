//determine whether fastq file is a 33 or 64 offset, store in qual_offset
//if scale is invalid, set qual_offset to -1 and return false
#include <cstdio>

bool fastq_file_offset(FILE * fastq_fh, int * qual_offset);
