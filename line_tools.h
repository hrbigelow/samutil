#ifndef _LINE_TOOLS_H
#define _LINE_TOOLS_H

#include <cstddef>
#include <cstdarg>

namespace LineTools
{
    void parse_fields(char const* format, 
                      char const* line,
                      char const* type_code,
                      void ** store);


    void parse_fields_aux(char const* format, 
                          char const* line,
                          char const* type_code,
                          void ** store,
                          ...);
    
    void parse_ordered_fields(char const* line,
                              char const* format, 
                              size_t const* order,
                              void ** store);
    
    int compare(void const* a, void const* b, char typecode);
    
    int compare_all(void ** af, void ** bf, 
                    char const* typecode, size_t num_fields);
};

#endif // _LINE_TOOLS_H
