//Convert Illumina GAPipeline 1.3.2 qseq format to FastQ format, as a filter
#include <algorithm>

int main(){

  int lane, tile, x, y, read;
  char machine[256];
  char sequence[1024];
  char qseq_string[1024];
  char fastq_string[1024];

  char qseq2fastq[256];
  for (int i=0; i < 256; ++i){
    qseq2fastq[i] = char(std::max(0,i-31));
  }

  char *q, *f;
  while (! feof(stdin)){
    fscanf(stdin, "%s %*i %i %i %i %i %*i %i %s %s %*i\n", machine, &lane, &tile, &x, &y, &read, sequence, qseq_string);
    for (q = qseq_string, f = fastq_string; *q != 0; ++q, ++f){ *f = qseq2fastq[int(*q)]; }
    *f = 0;
    fprintf(stdout, "@%s:%i:%i:%i:%i/%i\n%s\n+\n%s\n", machine, lane, tile, x, y, read, sequence, fastq_string);
  }
  
  return 0;
}
