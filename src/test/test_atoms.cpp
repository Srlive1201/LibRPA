#include "../atoms.h"
#include <iostream>

int main (int argc, char *argv[])
{
    using namespace std;
    double non_double;
    atom_mapping<double>::single_t a_double;
    atom_mapping<double>::pair_t_old apair_double;
    a_double[1] = 1.0;
    a_double[8] = -1.0;
    for (auto &a: get_atom(a_double))
        cout << a << endl;
    apair_double[1][1] = 1.0;
    apair_double[0][3] = 0.0;
    for (auto &a: get_atom_pair(apair_double))
        cout << a.first <<  " " << a.second << endl;

    return 0;
}
