#include <set>
#include <string>

struct samline_fragment_is_member
{
    std::set<std::string> id_list;
    samline_fragment_is_member(std::set<std::string> _id_list);
    bool operator()(char * samline);
};
