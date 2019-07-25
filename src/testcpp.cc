//clang 3.8.0

#include <iostream>

int main()
{
#ifdef MACRO
    std::cout << "Hello world!" << std::endl;
#else
    std::cerr << "Something wrong happened" << std::endl;
    return 1;
#endif
}
