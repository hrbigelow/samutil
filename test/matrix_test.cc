#include <cassert>
#include <cstdio>

#include "matrix_tools.h"

int main()
{
    int mat[] = {
        1, 0, 1, 2, 3,
        3, 2, 7, 1, 0,
        5, 7, 2, 5, 1
    };

    size_t dims[] = { 5, 3 };
    size_t ndims = 2;
    bool dims_to_keep[] = { true, false };

    int marg[5];

    MatrixTools::marginalize(mat, dims, ndims, dims_to_keep, marg);
    for (size_t i = 0; i != 5; ++i)
    {
        printf("%i\n", marg[i]);
    }

    int marg2[3];
    bool dims_to_keep2[] = { false, true };

    MatrixTools::marginalize(mat, dims, ndims, dims_to_keep2, marg2);
    for (size_t i = 0; i != 3; ++i)
    {
        printf("%i\n", marg2[i]);
    }

    //test the counters
    size_t dims3[] = { 5,8,2,10,8,4 };
    //size_t dims3[] = { 10, 10, 10, 10, 10, 10 };

    size_t counter[] = { 0, 0, 0, 0, 0, 0 };

    int keep[] = { 0, 1, 2, 3, 4, 5, -1 };
    
    size_t index_simple = 0;
    size_t index_cond;
    do
    {
        index_cond = MatrixTools::condense_digits(dims3, 6, counter, keep);
        assert(index_simple == index_cond);
        ++index_simple;

    } 
    while (MatrixTools::increment_counter(dims3, keep, counter));

}
