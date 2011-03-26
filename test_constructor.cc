#include <cstdlib>
#include <vector>

class Con
{
    int _a;
public:
    Con(int a)
    {
        Init(a);
    }
    Con(char const* as)
    {
        int a = atoi(as);
        Init(a);
    }
    void Init(int a)
    {
        this->_a = a;
    }
};


int main()
{
    std::vector<char **> bb;

    int a = 5;
    char * as = "5";
    Con c(a);
    Con d(as);
    return 0;
}
