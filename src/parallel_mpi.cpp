#include "parallel_mpi.h"

#include <algorithm>
#include <stdexcept>

namespace LIBRPA {

namespace MPI_Wrapper {

}  // namespace MPI_Wraper


} // namespace LIBRPA

Parallel_MPI::Parallel_MPI(void)
{
}

Parallel_MPI::~Parallel_MPI(void)
{
    // MPI_Finalize();
    // cout<<"Parallel_MPI Object is being deleted !"<<endl;
}

vector<double> Parallel_MPI::pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m)
{
    vector<double> pack;
    pack.push_back(-999);		const size_t index_size = pack.size()-1;
    for(const auto &I_pair:Cs_m)
    {
        size_t I=I_pair.first;
        for(const auto &J_pair:I_pair.second )
        {
            size_t J=J_pair.first;
            for(const auto &R_pair:J_pair.second)
            {
                const Vector3_Order<int> &R=R_pair.first;
                const matrix &m= *R_pair.second;
                pack.push_back(I);
                pack.push_back(J);
                pack.push_back(R.x);pack.push_back(R.y);pack.push_back(R.z);
                pack.push_back(m.nr);pack.push_back(m.nc);
                pack.insert(pack.end(), m.c, m.c+m.nr*m.nc);
            }
        }
    }
    pack[index_size] = pack.size() - (index_size+1);
    return pack;
}


// map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>>  Parallel_MPI::unpack_mat(vector<double> &pack)
// {
//     map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Cs;
//     auto ptr=pack.begin();
//     const size_t pack_size=ptr[0];
//     ptr+=1;
//     auto ptr_end=ptr+pack_size;
//     while(ptr<ptr_end)
//     {
//         const size_t I=ptr[0], J=ptr[1];
//         const Vector3_Order<int> R={ptr[2],ptr[3],ptr[4]};
//         const int nr=ptr[5], nc=ptr[6];
//         matrix &m=*Cs[I][J][R];
//         m.create(nr,nc);
//         copy(ptr+7,ptr+7+nr*nc,m.c);
//         ptr+=7+nr*nc;
//     }
//     return Cs;
// }
