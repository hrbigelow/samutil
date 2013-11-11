#include "gzip_tools.h"

#include <cstring>
#include <cstdlib>

#include <zlib.h>

// all from Mark Adler's pigz.c code
#define PUT2L(a,b) (*(a)=(b)&0xff,(a)[1]=(b)>>8)
#define PUT4L(a,b) (PUT2L(a,(b)&0xffff),PUT2L((a)+2,(b)>>16))
#define PUT4M(a,b) (*(a)=(b)>>24,(a)[1]=(b)>>16,(a)[2]=(b)>>8,(a)[3]=(b))


// fragment of put_header from pigz.c.  writes a gzip header
void put_gzip_header(char * name, time_t & mtime, 
                     int compression_level, FILE * out_fh)
{
    unsigned char head[10];

    head[0] = 31;
    head[1] = 139;
    head[2] = 8;                /* deflate */
    head[3] = name != NULL ? 8 : 0;
    PUT4L(head + 4, mtime);
    head[8] = compression_level == 1 ? 4 : 0;
    head[9] = 3;                /* unix */
    fwrite(head, 1, 10, out_fh);
    if (name != NULL)
    {
        fwrite(name, 1, strlen(name) + 1, out_fh);
    }
}



/* fragment of put_trailer from pigz.c.  writes a gzip trailer */
void put_gzip_trailer(unsigned long ulen, unsigned long check, FILE * out_fh)
{
    unsigned char tail[46];

    PUT4L(tail, check);
    PUT4L(tail + 4, ulen);
    fwrite(tail, 1, 8, out_fh);
}


unsigned long crc32_initial()
{
    return crc32(0L, Z_NULL, 0);
}


// return a crc32 checksum for this buffer chunk
unsigned long crc32_chunk(char const* buf, size_t buflen)
{
    unsigned long check = crc32(0L, Z_NULL, 0);
    size_t left = buflen;
    unsigned char const* next = reinterpret_cast<unsigned char const*>(buf);
    while (left > MAXP2) {
        check = crc32(check, next, MAXP2);
        left -= MAXP2;
        next += MAXP2;
    }
    check = crc32(check, next, (unsigned)left);
    return check;
}


// cut-and-paste from crc32.c
unsigned long gf2_matrix_times(unsigned long *mat, unsigned long vec)
{
    unsigned long sum;

    sum = 0;
    while (vec) {
        if (vec & 1)
            sum ^= *mat;
        vec >>= 1;
        mat++;
    }
    return sum;
}

// cut-and-paste from crc32.c
void gf2_matrix_square(unsigned long *square, unsigned long *mat)
{
    int n;

    for (n = 0; n < 32; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}


// cut and pasted from Mark Adler's pigz code.
unsigned long crc32_comb(unsigned long crc1, unsigned long crc2,
                         size_t len2)
{
    int n;
    unsigned long row;
    unsigned long even[32];     /* even-power-of-two zeros operator */
    unsigned long odd[32];      /* odd-power-of-two zeros operator */

    /* degenerate case */
    if (len2 == 0)
        return crc1;

    /* put operator for one zero bit in odd */
    odd[0] = 0xedb88320UL;          /* CRC-32 polynomial */
    row = 1;
    for (n = 1; n < 32; n++) {
        odd[n] = row;
        row <<= 1;
    }

    /* put operator for two zero bits in even */
    gf2_matrix_square(even, odd);

    /* put operator for four zero bits in odd */
    gf2_matrix_square(odd, even);

    /* apply len2 zeros to crc1 (first square will put the operator for one
       zero byte, eight zero bits, in even) */
    do {
        /* apply zeros operator for this bit of len2 */
        gf2_matrix_square(even, odd);
        if (len2 & 1)
            crc1 = gf2_matrix_times(even, crc1);
        len2 >>= 1;

        /* if no more bits set, then done */
        if (len2 == 0)
            break;

        /* another iteration of the loop with odd and even swapped */
        gf2_matrix_square(odd, even);
        if (len2 & 1)
            crc1 = gf2_matrix_times(odd, crc1);
        len2 >>= 1;

        /* if no more bits set, then done */
    } while (len2 != 0);

    /* return combined crc */
    crc1 ^= crc2;
    return crc1;
}
