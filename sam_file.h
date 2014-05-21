#ifndef _SAM_FILE_H
#define _SAM_FILE_H

#include <cstdio>

/*
  
 */

void SetToFirstDataLine(FILE ** sam_fh);

char * ReadAllocSAMHeader(FILE * sam_fh);

void PrintSAMHeader(FILE ** input_sam_fh, FILE * output_fh);

#endif // _SAM_FILE_H
