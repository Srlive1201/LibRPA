#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../core/coulmat.h"
#include "../core/exx.h"
#include "../io/global_io.h"
#include "../mpi/global_mpi.h"
#include "../utils/constants.h"

using namespace librpa_int;

static double exchange(const bool qstar, const bool clear_sector, const bool inversion = false)
{
    const auto& comm = global::mpi_comm_global_h;
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});
    pbc.set_irreducible_kgrids_kvec(3, 1, 1, {0, 0, 0, TWO_PI / 3, 0, 0},
                                    {{{0, 0, 0}}, {{1.0 / 3, 0, 0}, {-1.0 / 3, 0, 0}}});
    AtomicBasis basis(std::vector<std::size_t>{1});
    basis.set_l_shells({{0}});
    SymmetryContext ctx;
    ctx.set_crystal_structure(pbc.latvec, pbc.G, {{0, 0}}, {{0, {0, 0, 0}}});
    ctx.set_rspace_operations({SpaceGroupSymOp::IDENTITY});
    if (inversion)
    {
        SymmetryOperation inv;
        inv.rotation = Matrix3(-1, 0, 0, 0, -1, 0, 0, 0, -1);
        ctx.set_rspace_operations({SpaceGroupSymOp::IDENTITY, inv});
    }
    ctx.build_periodic_mappings(pbc, pbc.Rlist);
    ctx.build_rsh_rotations(
        {-1, 0, LIBRPA_ANGULAR_ORDER_ABS_PM, LIBRPA_RSH_COEFF_M_1, LIBRPA_RSH_COEFF_1_M}, 0);
    ctx.build_kstar_member_rotations(0);
    ctx.set_available();
    if (!inversion && ctx.count_irreducible_blocks() != pbc.Rlist.size())
        throw std::runtime_error("Identity fixture must not reduce real-space sectors");
    if (inversion && ctx.count_irreducible_blocks() >= pbc.Rlist.size())
        throw std::runtime_error("Inversion fixture must reduce real-space sectors");
    if (clear_sector)
    {
        ctx.irreducible_sector.clear();
        ctx.rspace_sector_stars.clear();
    }

    MeanField mf(1, 2, 1, 1);
    for (int ik = 0; ik < 2; ++ik)
    {
        mf.get_weight()[0](ik, 0) = 2.0 / 3.0;
        mf.get_eigenvals()[0](ik, 0) = -1.0;
        mf.get_eigenvectors()[0][0][ik].create(1, 1);
        mf.get_eigenvectors()[0][0][ik](0, 0) = 1.0;
    }
    Cs_LRI cs;
    atpair_k_cplx_mat_t vq;
    if (comm.myid == 0)
    {
        RI::Tensor<double> c({1, 1, 1});
        c(0, 0, 0) = 0.5;
        cs.data_libri[0][{0, {0, 0, 0}}] = c;
        for (const auto& q : pbc.klist)
        {
            auto v = std::make_shared<ComplexMatrix>(1, 1);
            (*v)(0, 0) = 2.0;
            vq[0][0][q] = v;
        }
    }
    bool replicated = false;
    const auto vr = FT_Vq(comm, basis, ctx, vq, pbc, true, qstar, &replicated);
    if (replicated != (qstar && !clear_sector))
        throw std::runtime_error("Unexpected Fourier output ownership receipt");
    // Two ranks form two 1x1 BLACS groups; a single rank keeps one 1x1 group.
    // The BLACS axis must stay at one process: a shape that spreads it leaves one
    // rank per k-point group with an empty local wave-function block, and the
    // symmetry-restore path of the current upstream rejects a stored wave
    // function whose size does not match such an empty descriptor block.
    const KPointBlacsProcessShape shape =
        comm.nprocs == 2 ? KPointBlacsProcessShape{2, 1} : KPointBlacsProcessShape{1, comm.nprocs};
    KPointBlacsParallelContext kctx(shape, comm.comm, 2);
    auto desc = kctx.create_array_desc(1, 1, 1, 1);
    Exx exx(mf, basis, pbc, ctx, kctx, kctx, desc, false, qstar);
    exx.build(LIBRPA_ROUTING_LIBRI, basis, cs, vr, replicated);
    double local = 0.0;
    for (const auto& s : exx.exx_IJR)
        for (const auto& b : s.second)
            for (const auto& k : b.second)
                for (const auto& i : k.second)
                    for (const auto& j : i.second)
                        for (const auto& r : j.second)
                            if (r.first == Vector3_Order<int>{0, 0, 0}) local += r.second(0, 0);
    double total = 0;
    comm.allreduce(&local, &total, 1, MPI_SUM);
    return total;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    global::init_global_mpi(MPI_COMM_WORLD);
    global::init_global_io(false, "stdout", false);
    int failure = 0;
    {
        const double direct = exchange(false, false);
        const double identity = exchange(true, false);
        const double full_sector = exchange(true, true);
        const double reduced = exchange(true, false, true);
        if (global::mpi_comm_global_h.is_root())
            std::cout << "OWNERSHIP_RESULT ranks=" << global::mpi_comm_global_h.nprocs
                      << " direct=" << direct << " identity=" << identity
                      << " full_sector=" << full_sector << " reduced=" << reduced << std::endl;
        failure = std::abs(direct) < 1e-8 || std::abs(identity - direct) > 1e-11 ||
                  std::abs(full_sector - direct) > 1e-11 || std::abs(reduced - direct) > 1e-11;
    }
    global::finalize_global_io();
    global::finalize_global_mpi();
    MPI_Finalize();
    return failure;
}
