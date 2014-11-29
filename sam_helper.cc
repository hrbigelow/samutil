#include "sam_helper.h"

#include <algorithm>

char base_to_complement[] =
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XTXGXXXCXXXXXXNXXXXXAXXXXXXXXXXX"
    "XtXgXXXcXXXXXXnXXXXXaXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";


// copies the reverse complement of a range defined by [begin, end)
// into rcomp.  rcomp must not overlap [begin, end) range
void reverse_comp(char const* begin, char const* end, char * rcomp)
{

    while (begin != end)
    {
        --end;
        *rcomp = base_to_complement[static_cast<int>(*end)];
        ++rcomp;
    }
}


// outputs the reverse complement of a range defined by [begin, end)
// into rcomp.  rcomp may be equal to begin or not, but it must contain
// enough space
void reverse_comp_inplace(char * begin, char * end)
{
    
    char * comp;
    for (comp = begin; comp != end; ++comp)
    {
        *comp = base_to_complement[static_cast<int>(*comp)];
    }
    std::reverse(begin, end);
}
