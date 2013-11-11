#include <time.h>
#include <cstdio>
#include <limits.h>

/* largest power of 2 that fits in an unsigned int -- used to limit requests
   to zlib functions that use unsigned int lengths */
#define MAXP2 (UINT_MAX - (UINT_MAX >> 1))


void put_gzip_header(char * name, time_t & mtime, 
                     int compression_level, FILE * out_fh);

void put_gzip_trailer(unsigned long ulen, unsigned long check, FILE * out_fh);


// return the running crc32 for the buffer provided, using the seed if
// this is not the first part of the buffer
unsigned long crc32_chunk(char const* buf, size_t buflen);


unsigned long gf2_matrix_times(unsigned long *mat, unsigned long vec);

void gf2_matrix_square(unsigned long *square, unsigned long *mat);

unsigned long crc32_comb(unsigned long crc1, unsigned long crc2,
                         size_t len2);

