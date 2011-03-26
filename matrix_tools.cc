#include "matrix_tools.h"

#include <algorithm>

namespace MatrixTools
{

    //recursive function to increment a counter, allowing masking certain digits
    bool increment_counter(size_t const* dim_sizes, 
                           int const* digit, //pointer to the digit that needs moving
                           size_t * counter)
    {
        counter[*digit]++;
        if (counter[*digit] == dim_sizes[*digit])
        {
            counter[*digit++] = 0;
            
            if (*digit == -1)
            {
                return false;
            }
            else
            {
                return 
                    increment_counter(dim_sizes, digit, counter);
            }
        }
        return true;
    }

    
    //recursive function to zero a counter, allowing to mask certain digits
    void zero_counter(size_t * counter, size_t num_digits, bool const* dims_to_zero)
    {
        for (size_t c = 0; c != num_digits; ++c)
        {
            if (dims_to_zero[c])
            {
                counter[c] = 0;
            }
        }
    }


    //calculate the numerical value of a mixed-base number,
    //represented as a set of digits of different size.  But, only use
    //<digits_to_keep> digits in the calculation For example, 18734,
    //with <digits_to_keep> = [0,2,3], would be 874, assuming
    //dim_sizes [10, 10, 10, 10, 10]
    size_t condense_digits(size_t const* dim_sizes,
                           size_t ndims,
                           size_t const* digits,
                           int *digits_to_keep)
    {
        int * last_digit = std::find(digits_to_keep, digits_to_keep + ndims, -1);
        size_t nkeep_dims = std::distance(digits_to_keep, last_digit);

        size_t num = digits[digits_to_keep[nkeep_dims - 1]];
        int di;
        for (size_t d = nkeep_dims - 1; d != 0; --d)
        {
            di = digits_to_keep[d - 1];
            num *= dim_sizes[di];
            num += digits[di];
        }
        return num;
    }


};
