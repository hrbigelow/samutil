#ifndef _MATRIX_TOOLS_INL_H
#define _MATRIX_TOOLS_INL_H

#include <algorithm>

namespace MatrixTools
{

    //marginalize out dimensions of a <num_dims> - dimensionl matrix called <source>
    //into <marginal>, keeping the dimensions <dims_to_keep>
    template <typename N>
    void marginalize(N const* source, 
                     size_t const* dim_sizes,
                     size_t ndims,
                     bool const* dims_to_keep,
                     N * marginal)
    {

        size_t * counter = new size_t[ndims];
        std::fill(counter, counter + ndims, 0);

        int * keep_digits = new int[ndims + 1];
        int * all_digits = new int[ndims + 1];

        size_t k = 0;
        size_t marg_size = 1;
        for (size_t d = 0; d != ndims; ++d)
        {
            if (dims_to_keep[d])
            {
                keep_digits[k++] = d;
                marg_size *= dim_sizes[d];
            }
            all_digits[d] = d;
        }

        std::fill(marginal, marginal + marg_size, 0);

        //sentinal value
        all_digits[ndims] = -1;
        keep_digits[k] = -1;

        size_t full_index = 0;
        do
        {
            size_t keep_index = condense_digits(dim_sizes, ndims, counter, keep_digits);
            marginal[keep_index] += source[full_index];
            ++full_index;
        }
        while (increment_counter(dim_sizes, all_digits, counter));
        
        delete counter;
        delete keep_digits;
        delete all_digits;

    }

}

#endif // _MATRIX_TOOLS_INL_H
