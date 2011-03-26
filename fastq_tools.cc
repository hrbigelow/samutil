//parse lines of a fastq file until all but one possible quality scale
//is eliminated
//return the offset

#include "fastq_tools.h"

//determine whether fastq file is a 33 or 64 offset, store in qual_offset
//if scale is invalid, set qual_offset to -1 and return false
bool fastq_file_offset(FILE * fastq_fh, int * qual_offset)
{
    bool file_types[] = { true, true, true, true };
    char const* min_qualcodes = "!;@B";
    char const* bound_qualcodes = "Jiii";

    int nfields_read;
    char quality_string[1024];
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
    while (! feof(fastq_fh)){
        nfields_read = 
            fscanf(fastq_fh, "%*s\n%*s\n%*s\n%s\n", quality_string);
        
        ++chunks;
        for (qual = quality_string; *qual != '\0'; ++qual)
        {
            for (size_t qc = 0; qc != 4; ++qc)
            {
                file_types[qc] = file_types[qc] 
                    && legal_values[qc][static_cast<size_t>(*qual)];
            }
        }
        if (file_types[0] &&
            ! (file_types[1] || file_types[2] || file_types[3]))
        {
            file_type = 0;
            *qual_offset = 33;
            break;
        }
        else if ((file_types[1] || file_types[2] || file_types[3]) &&
                 ! file_types[0])
        {
            file_type = 1;
            *qual_offset = 64;
            break;
        }
        if (chunks > 1000000)
        {
            break;
        }
            
    }
    return (file_type == 0 || file_type == 1);
}
