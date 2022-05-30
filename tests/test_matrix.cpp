#include "matrix.h"
#include <stdexcept>
#include <cstdlib>
#include <ctime>

int main (int argc, char *argv[])
{
    using namespace std;
    matrix m1(5, 6, true);
    srand(time(0));
    for (int i = 0; i < m1.size; i++)
        m1.c[i] = rand();

    matrix m2(m1);
    assert(m1.size == m2.size);
    for (int i = 0; i < m1.size; i++)
        if ( m1.c[i] != m2.c[i] ) throw logic_error("copy constructor fail");
    if (! (m1 == m2 )) throw logic_error("direct compare fail");
    return 0;

}

