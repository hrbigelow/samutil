#include "line_tools.h"

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <cstdlib>


/*
  how to get a scanf interface that guesses the number of
  fields
 */

void LineTools::parse_fields(char const* format, 
                             char const* line,
                             char const* type_code,
                             void ** store)
{
    void ** s = store;

    switch(strlen(type_code))
    {
    case 0: return;
        break;
        
    case 1: return LineTools::parse_fields_aux(format, line, type_code, s, s[0]);
        break;

    case 2: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1]);
        break;

    case 3: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2]);
        break;

    case 4: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3]);
        break;

    case 5: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3], s[4]);
        break;

    case 6: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3], s[4], s[5]);
        break;

    case 7: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3], s[4], s[5], 
                                               s[6]);
        break;

    case 8: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3], s[4], s[5], 
                                               s[6], s[7]);
        break;

    case 9: return LineTools::parse_fields_aux(format, line, type_code, store,
                                               s[0], s[1], s[2], s[3], s[4], s[5], 
                                               s[6], s[7], s[8]);
        break;

    case 10: return LineTools::parse_fields_aux(format, line, type_code, store,
                                                s[0], s[1], s[2], s[3], s[4], s[5], 
                                                s[6], s[7], s[8], s[9]);
        break;

    default:
        fprintf(stderr, "Error: LineTools::parse_fields: "
                "Cannot handle more than 10 fields\n");
        exit(1);
        break;
    }
}


//parse and load fields from an input line according to a format
//string convert string values to int, string, or float.
//do this generically w.r.t the number and formats of the fields.
//the output, <store> is intended to be used generically with the
//matching format specifiers
void LineTools::parse_fields_aux(char const* format, 
                                 char const* line,
                                 char const* type_code,
                                 void ** store,
                                 ...)
{

    va_list args;
    va_start(args, store);
    size_t count = static_cast<size_t>(vsscanf(line, format, args));

    for (size_t t = 0; t != count; ++t)
    {
        switch(type_code[t])
        {
        case 'i':
            memcpy(store[t], va_arg(args, int const*), sizeof(int));
            break;
        case 'f':
            memcpy(store[t], va_arg(args, float const*), sizeof(float));
            break;
        case 's':
            {
                char const* arg = va_arg(args, char const*); 
                memcpy(store[t], arg, strlen(arg));
                break;
            }
        default:
            fprintf(stderr, "Error: unrecognized typecode '%c'\n", type_code[t]);
            exit(1);
            break;
        }
    }
    va_end(args);
}


//parse a set of ordered fields, loading them into store.  'order'.
//Assume the number N of active format fields (not with %*s, for
//example) will fit in 'store', and order[0..N] contains numbers 0 through N
//stores 
/*
void LineTools::parse_ordered_fields(char const* line,
                                     char const* format, 
                                     size_t const* order,
                                     void ** store)
{
    va_list args;
    va_start(args, format);
    int count = vsscanf(line, format, args);
    va_end(args);
    size_t a = 0;
    char const* arg;
    while ((arg = va_arg(args, char const*)) != NULL)
    {
        if (order[a] >= count)
        {
            fprintf(stderr, 
                    "error: LineTools::parse_ordered_fields: "
                    "ordering array contains"
                    "out of bounds element %Zu at %Zu\n", order[a], a);
            exit(1);
        }
        strcpy(store[order[a]], arg);
        ++a;
    }
}
*/


int LineTools::compare(void const* a, void const* b, char typecode)
{
    switch (typecode)
    {
    case 'i':
        {
            int anum = *static_cast<int const*>(a);
            int bnum = *static_cast<int const*>(b);
            return anum < bnum ? -1 : (anum == bnum ? 0 : 1);
            break;
        }
    case 'f':
        {
            float anum = *static_cast<float const*>(a);
            float bnum = *static_cast<float const*>(b);
            return anum < bnum ? -1 : (anum == bnum ? 0 : 1);
            break;
        }
    case 's':
        {
            char const* astr = static_cast<char const*>(a);
            char const* bstr = static_cast<char const*>(b);
            return strcmp(astr, bstr);
            break;
        }
    default:
        fprintf(stderr, "Error: unrecognized typecode '%c'\n", typecode);
        exit(1);
        break;
    }
}


int LineTools::compare_all(void ** af, void ** bf, 
                           char const* typecode, size_t num_fields)
{
    int cmp = 0;
    for (size_t f = 0; f != num_fields; ++f)
    {
        cmp = compare(af[f], bf[f], typecode[f]);
        if (cmp != 0)
        {
            break;
        }
    }
    return cmp;
}
