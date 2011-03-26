#include <set>

struct switch_less
{
    bool switchy;
    switch_less(bool _switchy = true) : switchy(_switchy) { }
    //switch_less(switch_less const& sl) : switchy(sl.switchy) { }

    bool operator()(int a, int b)
    {
        if (switchy)
        {
            return a < b;
        }
        else
        {
            return b < a;
        }
    }
};


int main(int argc, char ** argv)
{

    switch_less sl(true);

    std::set<int, switch_less> switch_set(sl);
    switch_set.insert(3);
    switch_set.insert(5);
    for (std::set<int, switch_less>::iterator it = switch_set.begin();
         it != switch_set.end(); ++it)
    {
        fprintf(stdout, "%i\n", *it);
    }
    return 0;
};
