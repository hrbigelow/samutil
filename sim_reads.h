#ifndef _SIM_READS_H
#define _SIM_READS_H

#include <cstdlib>

int sim_reads_usage(char const* ddef, size_t ldef, 
                    char const* sdef,
                    char cdef, char udef,
                    char const* mdef, size_t rdef);

int main_sim_reads(int argc, char ** argv);

#endif // _SIM_READS_H
