//parse lines of a fastq file until all but one possible quality scale
//is eliminated
//return the offset

#include "fastq_tools.h"
#include <algorithm>
#include <cstring>

/*
//determine whether fastq file is a 33 or 64 offset, store in qual_offset
//if scale is invalid, set qual_offset to -1 and return false
bool fastq_file_offset(FILE * fastq_fh, int * qual_offset,
                       char ** file_type_strings)
{
    bool file_types[] = { true, true, true, true };
    char const* min_qualcodes = "!;@B";
    char const* bound_qualcodes = "Jiii";

    int nfields_read;
    char qualstr[1024];
    char * qual;

    bool legal_values[4][256];
    for (size_t st = 0; st != 4; ++st)
    {
        for (size_t vt = 0; vt != 256; ++vt)
        {
            legal_values[st][vt] = 
                min_qualcodes[st] <= static_cast<char>(vt) 
                && static_cast<char>(vt) < bound_qualcodes[st];
        }
    }
    *qual_offset = -1;
    size_t file_type = 5;
    size_t chunks = 0;

    char maxc = 0, minc = 254;
    char s1, s2, s3;

    while (! feof(fastq_fh)){
        nfields_read = 
            fscanf(fastq_fh, "%*s\n%*s\n%*s\n%s\n", qualstr);

        size_t s = strlen(qualstr);
        s1 = qualstr[0];
        s2 = qualstr[(s-1)/2];
        s3 = qualstr[s-1];

        // sample first, last, and middle characters
        maxc = std::max(std::max(s1, s2), std::max(s3, maxc));
        minc = std::min(std::min(s1, s2), std::min(s3, minc));

        ++chunks;
        // for (qual = qualstr; *qual != '\0'; ++qual)
        // {
        //     for (size_t qc = 0; qc != 4; ++qc)
        //     {
        //         file_types[qc] = file_types[qc] 
        //             && legal_values[qc][static_cast<size_t>(*qual)];
        //     }
        // }
    }
    if (file_types[0] && ! file_types[1])
    {
        file_type = 0;
        *qual_offset = 33;
    }
    else if (! file_types[0] && file_types[1] && file_types[2] && file_types[3])
    {
        file_type = 1;
        *qual_offset = 64;
    }

    char const* file_type_names = "SXIJ";
    for (size_t i = 0; i != 4; ++i)
    {
        file_types[i] =
            legal_values[qc][static_cast<size_t>(maxc)]
            && legal_values[qc][static_cast<size_t>(minc)];

        (*file_type_strings)[i] = file_types[i] ? file_type_names[i] : '-';
    }

    return (file_type == 0 || file_type == 1);
}
*/


void fastq_extreme_chars(FILE * fastq_fh, char * minc, char * maxc)
{
    int nfields_read;
    char qualstr[1024];

    *maxc = static_cast<char>(0);
    *minc = static_cast<char>(127);
    char s1, s2, s3;

    while (! feof(fastq_fh)){
        nfields_read = 
            fscanf(fastq_fh, "%*[^\n]\n%*s\n%*s\n%s\n", qualstr);

        size_t s = strlen(qualstr);
        s1 = qualstr[0];
        s2 = qualstr[(s-1)/2];
        s3 = qualstr[s-1];

        // sample first, last, and middle characters
        *maxc = std::max(std::max(s1, s2), std::max(s3, *maxc));
        *minc = std::min(std::min(s1, s2), std::min(s3, *minc));
        
    }
}
