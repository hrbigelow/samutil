#include <gmp.h>

void int_to_Hilbert(mpz_t const& index, size_t nchunks, uint64_t * coords, size_t ndim);
void Hilbert_to_int(uint64_t const* coords, size_t ndim, size_t nchunks, mpz_t & index);

void unpack_index(mpz_t const& index, size_t ndim, int * chunks, size_t nchunks);
void pack_index(int const* chunks, size_t nchunks, size_t ndim, mpz_t & index);
void initial_start_end(size_t nchunks, size_t ndim, int * start, int * end);
void unpack_coords(uint64_t const* coords, size_t ndim, int * unpacked, size_t nchunks);
void pack_coords(int const* chunks, size_t nchunks, uint64_t * packed, size_t ndim);
int gray_encode(int const& bn);
int gray_decode(int const& gray);
int gray_encode_travel(int start, int end, int mask, int i);
int gray_decode_travel(int start, int end, int mask, int gray);
void child_start_end(int parent_start, int parent_end, int mask, int i,
                     int * child_start, int * child_end);
