#ifndef PARALLEL_MPI_H
#define PARALLEL_MPI_H

#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
#include <mpi.h>
#include <string>
#include <cassert>
#include <set>
#include <vector>
#include <utility>
#include <map>
#include <memory>

using std::vector;
using std::map;
using std::pair;

class Parallel_MPI
{
public:
    static vector<double> pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m);
    // static map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> unpack_mat(vector<double> &pack);

    Parallel_MPI();
    ~Parallel_MPI();

    // void set_blacs_parameters();
    // void set_blacs_mat(
    //     int *desc, int &loc_row, int &loc_col, 
    //     const int tot_row, const int tot_col,
    //     const int row_blk=1, const int col_blk=1 );
};

#endif
