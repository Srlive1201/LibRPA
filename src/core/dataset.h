#pragma once
#include "../mpi/base_mpi.h"

namespace librpa_int
{

#define is_null_dataset_ptr(p) p == nullptr

class Dataset
{
public:
    MpiCommHandler comm_h_;

    Dataset(MPI_Comm comm): comm_h_(comm) {}

    void free() {};
    ~Dataset() { free(); }
};

}
