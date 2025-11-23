#pragma once
#include "interface/mpi.h"

namespace librpa_int
{

#define is_null_dataset_ptr(p) p == nullptr

class Dataset
{
public:
    MPI_Comm comm_;

    Dataset(MPI_Comm comm): comm_(comm) {}

    void free() {};
    ~Dataset() { free(); }
};

}
