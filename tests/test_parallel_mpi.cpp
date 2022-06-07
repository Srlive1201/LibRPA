#include <mpi.h>
#include <stdexcept>
#include <stdexcept>
#include "parallel_mpi.h"

int main (int argc, char *argv[])
{
    para_mpi.mpi_init(argc, argv);
    if ( para_mpi.get_size() != 4 )
        throw invalid_argument("test imposes 4 MPI processes");

    const int myid = para_mpi.get_myid();
    const int size = para_mpi.get_size();

    // single index dispatcher
    // balanced sequential case
    //   proc 0: 0, 1, 2
    //   proc 1: 3, 4, 5
    //   proc 2: 6, 7, 8
    //   proc 3: 9, 10, 11
    auto ilist = dispatcher(0, 12, myid, size, true);
    /* for (auto &id: ilist) */
    /*     printf("myid = %d, id = %u\n", myid, id); */
    assert(ilist.size() == 3);
    for ( int i = 0; i != ilist.size(); i++ )
        assert( ilist[i] == myid * 3 + i);

    // balanced non-sequential case
    //   proc 0: 0, 4, 8
    //   proc 1: 1, 5, 9
    //   proc 2: 2, 6, 10
    //   proc 3: 3, 7, 11
    ilist = dispatcher(0, 12, myid, size, false);
    assert (ilist.size() == 3);
    for ( int i = 0; i != ilist.size(); i++ )
        assert( ilist[i] == myid + i * 4);

    // inbalanced non-sequential dispatching
    //   proc 0: 0, 4
    //   proc 1: 1
    //   proc 2: 2
    //   proc 3: 3
    ilist = dispatcher(0, 5, myid, size, false);
    /* for (auto &id: ilist) */
    /*     printf("myid = %d, id = %u\n", myid, id); */
    if (myid == 0)
    {
        assert(ilist.size() == 2);
        assert(ilist[0] == 0);
        assert(ilist[1] == 4);
    }
    else
    {
        assert(ilist.size() == 1);
        assert(ilist[0] == myid);
    }

    // inbalanced sequential dispatching
    //   proc 0: 0, 1
    //   proc 1: 2, 3
    //   proc 2: 4
    //   proc 3: 5
    ilist = dispatcher(0, 6, myid, size, true);
    /* for (auto &id: ilist) */
    /*     printf("myid = %d, id = %u\n", myid, id); */
    if ((myid == 0) || (myid == 1))
    {
        assert(ilist.size() == 2);
        assert(ilist[0] == myid * 2);
        assert(ilist[1] == myid * 2 + 1);
    }
    else
    {
        assert(ilist.size() == 1);
        assert(ilist[0] == myid + 2);
    }

    // =====================================
    // double index dispatcher
    // balanced sequential, 1st faster
    //   proc 0: (0, 1), (1, 1)
    //   proc 1: (0, 2), (1, 2)
    //   proc 2: (0, 3), (1, 3)
    //   proc 3: (0, 4), (1, 4)
    auto i2list = dispatcher(0, 2, 1, 5, myid, size, true, true);
    assert(i2list.size() == 2);
    for ( int i = 0; i != 2; i++)
    {
        assert(i2list[i].first == i);
        assert(i2list[i].second == myid + 1);
    }

    // balanced sequential, 2nd faster
    //   proc 0: (0, 1), (0, 2)
    //   proc 1: (0, 3), (0, 4)
    //   proc 2: (1, 1), (1, 2)
    //   proc 3: (1, 3), (1, 4)
    i2list = dispatcher(0, 2, 1, 5, myid, size, true, false);
    assert(i2list.size() == 2);
    for ( int i = 0; i != 2; i++)
    {
        assert(i2list[i].first == myid / 2);
        assert(i2list[i].second == (myid % 2) * 2 + i + 1);
    }
    return 0;
    // balanced non-sequential, 1st faster
    //   proc 0: (0, 1), (0, 3)
    //   proc 1: (1, 1), (1, 3)
    //   proc 2: (0, 2), (0, 4)
    //   proc 3: (1, 2), (1, 4)
    i2list = dispatcher(0, 2, 1, 5, myid, size, false, true);
    assert(i2list.size() == 2);
    for ( int i = 0; i != 2; i++)
    {
        assert(i2list[i].first == myid % 2);
        assert(i2list[i].second == 1 + (myid / 2) + i * 2);
    }
    // balanced non-sequential, 2nd faster
    //   proc 0: (0, 1), (2, 1)
    //   proc 1: (0, 2), (2, 2)
    //   proc 2: (1, 1), (3, 1)
    //   proc 3: (1, 2), (3, 2)
    i2list = dispatcher(0, 4, 1, 3, myid, size, false, false);
    assert(i2list.size() == 2);
    for ( int i = 0; i != 2; i++)
    {
        assert(i2list[i].first == myid / 2 + i * 2);
        assert(i2list[i].second == myid % 2 + 1);
    }
    return 0;

    // inbalanced sequential, 1st faster
    i2list = dispatcher(0, 3, 1, 3, myid, size, true, true);
    //   proc 0: (0, 1), (1, 1)
    //   proc 1: (2, 1), (0, 2)
    //   proc 2: (1, 2)
    //   proc 3: (2, 2)
    if ( (myid == 0) || (myid == 1) )
    {
        assert(i2list.size() == 2);
        assert(i2list[0].first == myid * 2);
        assert(i2list[0].second == 1);
        assert(i2list[1].first == 1 - myid);
        assert(i2list[1].second == myid + 1);
    }
    else
    {
        assert(i2list.size() == 1);
        assert(i2list[0].first == myid - 1);
        assert(i2list[0].second == 2);
    }
    // inbalanced sequential, 2nd faster
    i2list = dispatcher(0, 2, 2, 5, myid, size, true, false);
    //   proc 0: (0, 2), (0, 3)
    //   proc 1: (0, 4), (1, 2)
    //   proc 2: (1, 3)
    //   proc 3: (1, 4)
    if ( (myid == 0) || (myid == 1) )
    {
        assert(i2list.size() == 2);
        assert(i2list[0].first == 0);
        assert(i2list[0].second == (myid+1)*2);
        assert(i2list[1].first == myid);
        assert(i2list[1].second == 3 - myid);
    }
    else
    {
        assert(i2list.size() == 1);
        assert(i2list[0].first == 1);
        assert(i2list[0].second == myid + 1);
    }

    // inbalanced non-sequential, 1st faster
    i2list = dispatcher(-1, 2, 1, 3, myid, size, false, true);
    //   proc 0: (-1, 1), (1, 1)
    //   proc 1: (-1, 2), (1, 2)
    //   proc 2: (0, 1)
    //   proc 3: (0, 2)
    if ( (myid == 0) || (myid == 1) )
    {
        assert(i2list.size() == 2);
        assert(i2list[0].first == -1);
        assert(i2list[0].second == myid + 1);
        assert(i2list[1].first == 1);
        assert(i2list[1].second == myid + 1);
    }
    else
    {
        assert(i2list.size() == 1);
        assert(i2list[0].first == 0);
        assert(i2list[0].second == myid - 1);
    }
    // inbalanced non-sequential, 2nd faster
    i2list = dispatcher(0, 3, -2, 0, myid, size, true, false);
    //   proc 0: (0, -2), (2, -2)
    //   proc 1: (0, -1), (2, -1)
    //   proc 2: (1, -2)
    //   proc 3: (1, -1)
    if ( (myid == 0) || (myid == 1) )
    {
        assert(i2list.size() == 2);
        assert(i2list[0].first == 0);
        assert(i2list[0].second == myid-2);
        assert(i2list[1].first == 2);
        assert(i2list[1].second == myid-2);
    }
    else
    {
        assert(i2list.size() == 1);
        assert(i2list[0].first == 1);
        assert(i2list[0].second == myid - 4);
    }
}
