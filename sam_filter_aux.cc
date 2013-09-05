#include "sam_filter_aux.h"

#include <set>
#include <string>
#include <string.h>

//unary function for computing index
samline_fragment_is_member::samline_fragment_is_member(std::set<std::string> _id_list) :
    id_list(_id_list) { }

bool samline_fragment_is_member::operator()(char * samline)
{
    char id[256];
    char * lptr = samline;
    lptr = strchr(samline, '\t');
    strncpy(id, samline, lptr - samline);
    id[lptr - samline] = '\0';
    // printf("%s\n", id);
    return this->id_list.find(id) != this->id_list.end();
}
