#ifndef _MATRIX_TOOLS_H
#define _MATRIX_TOOLS_H

#include <cstddef> // for size_t, ha!

namespace MatrixTools
{
    bool increment_counter(size_t const* dim_sizes, 
                           int const* digit, //pointer to the digit that needs moving
                           size_t * counter);
    
    void zero_counter(size_t * counter, size_t num_digits, bool const* dims_to_zero);
    
    size_t condense_digits(size_t const* dim_sizes,
                           size_t ndims,
                           size_t const* digits,
                           int * digits_to_keep);
    
};

#include "matrix_tools-inl.h"

#endif // _MATRIX_TOOLS_H
