#include "time_tools.h"

// miscellaneous helper function
size_t elapsed_ms(timespec & beg, timespec & end)
{
    return
        ((end.tv_sec * 1000000000 + end.tv_nsec) - 
         (beg.tv_sec * 1000000000 + beg.tv_nsec)) / 1000000;
}
