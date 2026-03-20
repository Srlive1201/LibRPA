// API headers
#include "librpa_enums.h"
#include "librpa_func.h"
#include "librpa_input.h"

// Standard headers
#include <cstring>
#include <ios>
#include <memory>
#include <vector>

// Internal headers
#include "../io/global_io.h"
#include "../io/stl_io_helper.h"
#include "../math/matrix.h"
#include "../math/vector3_order.h"
#include "../utils/error.h"
#include "instance_manager.h"

// External headers and stubs
#ifdef LIBRPA_USE_LIBRI
#include <initializer_list>
#include <RI/global/Tensor.h>
#else
#include "../utils/libri_stub.h"
#endif

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_scf_dimension,
                     int nspins, int nkpts, int nstates, int nbasis)
{
    using librpa_int::api::get_dataset_instance;
    using librpa_int::global::mpi_comm_global_h;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);

    auto &meanfield = pds->mf;
    int st_ib = 0;
    int nb_local = nstates;
    int st_iao = 0;
    int nao_local = nbasis;
    meanfield.set(nspins, nkpts, nstates, nbasis,
                  st_ib, nb_local, st_iao, nao_local);
    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        lib_printf("Mean-field dimensions set:\n");
        lib_printf("| number of spins    : %d\n", meanfield.get_n_spins());
        lib_printf("| number of k-points : %d\n", meanfield.get_n_kpoints());
        lib_printf("| number of bands    : %d\n", meanfield.get_n_bands());
        lib_printf("| number of NAOs     : %d\n", meanfield.get_n_aos());
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_wg_ekb_efermi, int nspins, int nkpts, int nstates, const double* wg,
                     const double* ekb, double efermi)
{
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &meanfield = pds->mf;

    meanfield.get_efermi() = efermi;
    auto& eskb = meanfield.get_eigenvals();
    auto& swg = meanfield.get_weight();
    int length_kb = nkpts * nstates;
    for (int is = 0; is != nspins; is++)
    {
        memcpy(eskb[is].c, ekb + length_kb * is, length_kb * sizeof(double));
        memcpy(swg[is].c, wg + length_kb * is, length_kb * sizeof(double));
        // wg[is](k_index, i) = stod(ws) / n_kpoints; // different with abacus!
        swg[is] *= (1.0 / nkpts);
    }

    double emin, emax;
    pds->mf.get_E_min_max(emin, emax);

    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        lib_printf("Mean-field eigenvalues and occupation numbers set:\n");
        lib_printf("| Minimal transition energy (Ha): %f\n", emin);
        lib_printf("| Maximal transition energy (Ha): %f\n", emax);
        lib_printf("| Fermi level               (Ha): %f\n", efermi);
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_wfc, int ispin, int ik,
                     int nstates_local, int nbasis_local,
                     const double *wfc_real, const double *wfc_imag)
{
    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &meanfield = pds->mf;

    auto& wfc = meanfield.get_eigenvectors()[ispin][ik];
    wfc.create(nstates_local, nbasis_local);
    const size_t n = meanfield.get_n_bands() * meanfield.get_n_aos();
    for (size_t i = 0; i < n; i++)
    {
        wfc.c[i] = std::complex<double>(wfc_real[i], wfc_imag[i]);
    }
    // std::cout << "Maxabs: " << wfc.get_max_abs() << std::endl;
    librpa_int::global::ofs_myid
        << "Wave-function set : ispin = " << ispin << " ik = " << ik
        << " nstates_local = " << nstates_local << " nbasis_local = " << nbasis_local
        << std::endl;
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ao_basis_wfc, const int natoms, const size_t *nbs_wfc)
{
    using librpa_int::global::lib_printf;

    std::vector<size_t> nbs(natoms);
    for (int i = 0; i < natoms; i++) nbs[i] = librpa_int::as_size(nbs_wfc[i]);

    auto pds = librpa_int::api::get_dataset_instance(h);
    pds->basis_wfc.set(nbs);
    pds->desc_wfc.reset_handler(pds->blacs_h);
    const auto n = pds->basis_wfc.nb_total;
    pds->desc_wfc.init_1b1p(n, n, 0, 0);

    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        lib_printf("Wave-function basis functions set:\n");
        lib_printf("| total number of basis: %lu\n", pds->basis_wfc.nb_total);
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ao_basis_aux, int natoms, const size_t *nbs_aux)
{
    using librpa_int::global::lib_printf;

    std::vector<size_t> nbs(natoms);
    librpa_int::global::ofs_myid << "Parsing basis: ";
    for (int i = 0; i < natoms; i++)
    {
        nbs[i] = librpa_int::as_size(nbs_aux[i]);
        librpa_int::global::ofs_myid << nbs[i] << " ";
    }
    librpa_int::global::ofs_myid << std::endl;

    auto pds = librpa_int::api::get_dataset_instance(h);
    pds->basis_aux.set(nbs);

    // After auxiliary basis is set, we can initialize the global (continous) array descriptor for N_abf size basis.
    pds->desc_abf.reset_handler(pds->blacs_h);
    const auto n = pds->basis_aux.nb_total;
    pds->desc_abf.init_1b1p(n, n, 0, 0);

    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        lib_printf("Auxiliary basis functions set:\n");
        lib_printf("| total number of basis: %lu\n", pds->basis_aux.nb_total);
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_latvec_and_G, const double lat_mat[9], const double G_mat[9])
{
    using std::cout;
    using std::endl;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &pbc = pds->pbc;

    std::vector<double> latt(lat_mat, lat_mat + 9);
    std::vector<double> recp(G_mat, G_mat + 9);

    pbc.set_latvec_and_G(latt, recp);

    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        cout << "Lattice vectors (Bohr): latt" << endl;
        pbc.latvec.print(16);
        cout << "Reciprocal lattice vectors (2PI Bohr^-1): G" << endl;
        pbc.G.print(16);
        const auto iden_test = pbc.latvec * pbc.G.Transpose();
        cout << "Consistency check: latt * G.T" << endl;
        iden_test.print(16);
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_atoms, int natoms, const int *types, const double *posi_cart)
{
    using std::cout;
    using std::endl;
    using librpa_int::coord_t;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &pbc = pds->pbc;
    auto &atoms = pds->atoms;

    std::vector<int> v_types(types, types + natoms);
    std::vector<coord_t> v_coords(natoms);
    for (int i = 0; i < natoms; i++)
    {
        v_coords[i][0] = posi_cart[3 * i];
        v_coords[i][1] = posi_cart[3 * i + 1];
        v_coords[i][2] = posi_cart[3 * i + 2];
    }

    // Set the fractional coordinates as well if the lattice has been set manually
    if (pbc.is_latt_set())
    {
        atoms.set(v_types, v_coords, pbc.latvec);
        const auto &coords = atoms.coords;
        const auto &coords_frac = atoms.coords_frac;
        pds->comm_h.barrier();
        if (pds->comm_h.is_root())
        {
            cout << "Atom positions read (Cartesian in Bohr | fractional):" << endl;
            for (int i = 0; i != natoms; i++)
            {
                const auto i_at = librpa_int::as_atom(i);
                const auto &c = coords.at(i_at);
                const auto &cf = coords_frac.at(i_at);
                lib_printf("ia %4d: %12.7f %12.7f %12.7f | %12.7f %12.7f %12.7f\n", i + 1, c[0], c[1],
                        c[2], cf[0], cf[1], cf[2]);
            }
        }
        pds->comm_h.barrier();
    }
    else
    {
        atoms.set(v_types, v_coords);
        const auto &coords = atoms.coords;
        pds->comm_h.barrier();
        if (pds->comm_h.is_root())
        {
            cout << "Atom positions read (Cartesian in Bohr, fractional not set due to uninitialized lattice):" << endl;
            for (int i = 0; i != natoms; i++)
            {
                const auto i_at = librpa_int::as_atom(i);
                const auto &c = coords.at(i_at);
                lib_printf("ia %4d: %12.7f %12.7f %12.7f\n", i + 1, c[0], c[1], c[2]);
            }
        }
        pds->comm_h.barrier();
    }

}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_kgrids_kvec, int nk1, int nk2, int nk3, const double* kvecs)
{
    using librpa_int::global::lib_printf;
    using std::cout;
    using std::endl;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &pbc = pds->pbc;

    const int nkpts = nk1 * nk2 * nk3;
    std::vector<double> v_kvecs(kvecs, kvecs + 3 * nkpts);

    pbc.set_kgrids_kvec(nk1, nk2, nk3, v_kvecs);

    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        librpa_int::global::lib_printf("kgrids: %3d %3d %3d\n", pbc.period.x, pbc.period.y, pbc.period.z);
        cout << "k-points read (Cartesian in 2Pi Bohr^-1 | fractional):" << endl;

        const auto &klist = pbc.klist;
        const auto &kfrac_list = pbc.kfrac_list;
        for (int ik = 0; ik != nkpts; ik++)
        {
            lib_printf("ik %4d: %12.7f %12.7f %12.7f | %12.7f %12.7f %12.7f\n",
                    ik+1, klist[ik].x, klist[ik].y, klist[ik].z,
                    kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
        }

        const auto &Rlist = pbc.Rlist;
        cout << "R-points to compute:" << endl;
        for (int iR = 0; iR != nkpts; iR++)
        {
            lib_printf("%4d: %3d %3d %3d\n", iR+1, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
        }
        cout << endl;
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ibz_mapping, int nkpts, const int* map_ibzk)
{
    using std::cout;
    using std::endl;
    using namespace librpa_int;  // for STL io
    using namespace librpa_int::global;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &pbc = pds->pbc;
    assert(nkpts == pbc.get_n_cells_bvk());

    std::vector<int> map(map_ibzk, map_ibzk + nkpts);
    pbc.set_ibz_mapping(map, {});
    // ofs_myid << map << std::endl;
    pds->comm_h.barrier();
    if (pds->comm_h.is_root())
    {
        const int nkpt = pbc.irk_point_id_mapping.size();
        cout << "Irreducible Brillouin mapping:" << endl;
        lib_printf("%4s: %12s %12s %12s    %4s\n",
                   "ik", "k_1", "k_2", "k_3", "ibzk");
        for (int ik = 0; ik < nkpt; ik++)
        {
            int ibzk = pbc.irk_point_id_mapping[ik];
            if (ibzk == ik)
            {
                lib_printf("%4d: %12.7f %12.7f %12.7f\n",
                           ik + 1, pbc.kfrac_list[ik].x, pbc.kfrac_list[ik].y, pbc.kfrac_list[ik].z);
            }
            else
            {
                lib_printf("%4d: %12.7f %12.7f %12.7f -> %4d\n",
                           ik + 1, pbc.kfrac_list[ik].x, pbc.kfrac_list[ik].y, pbc.kfrac_list[ik].z, ibzk + 1);
            }
        }
    }
    pds->comm_h.barrier();
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_lri_coeff, LibrpaParallelRouting routing, int I, int J,
                     int nbasis_i, int nbasis_j, int naux_mu, const int R[3], const double* Cs_in)
{
    using std::endl;
    using namespace librpa_int;
    using librpa_int::matrix;
    using librpa_int::Vector3_Order;
    using librpa_int::as_size;
    using librpa_int::global::ofs_myid;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &cs_data = pds->cs_data;
    if (!pds->basis_wfc.initialized())
        throw LIBRPA_RUNTIME_ERROR("wave function basis not set, call (librpa_)set_ao_basis_wfc first");
    if (!pds->basis_aux.initialized())
        throw LIBRPA_RUNTIME_ERROR("auxiliary basis not set, call (librpa_)set_ao_basis_aux first");

    const size_t cs_size = nbasis_i * nbasis_j * naux_mu;
    const size_t n_ij = nbasis_i * nbasis_j;

    // ofs_myid << "Parsing I J " << I << " " << J << " "
    //          << "ds->basis_aux " << pds->basis_aux.nbs_ << " "
    //          << "ds->basis_wfc " << pds->basis_wfc.nbs_ << endl;
    // ofs_myid << "ds->basis_aux[I] == naux_mu ? "
    //          << pds->basis_aux[I] << " " << as_size(naux_mu) << " "
    //          << std::boolalpha << (pds->basis_aux[I] == as_size(naux_mu))
    //          << endl;
    // ofs_myid << "ds->basis_wfc[I] == nbasis_i ? "
    //          << pds->basis_wfc[I] << " " << as_size(nbasis_i) << " "
    //          << std::boolalpha << (pds->basis_wfc[I] == as_size(nbasis_i))
    //          << endl;
    // ofs_myid << "ds->basis_wfc[J] == nbasis_j ? "
    //          << pds->basis_wfc[J] << " " << as_size(nbasis_j) << " "
    //          << std::boolalpha << (pds->basis_wfc[J] == as_size(nbasis_j))
    //          << endl;

    assert(pds->basis_aux[I] == as_size(naux_mu));
    assert(pds->basis_wfc[I] == as_size(nbasis_i));
    assert(pds->basis_wfc[J] == as_size(nbasis_j));

    if (routing == LibrpaParallelRouting::LIBRI)
    {
        const std::array<int, 3> Ra{R[0], R[1], R[2]};
        // RI tensor uses ABF as slowest index, so we need transpose the input data.
        auto data = std::make_shared<std::valarray<double>>(cs_size);
        // Unless n_cols >> n_rows, cache-friendly reading (by row in C/C++) is more efficient
        for (size_t i_row = 0; i_row < n_ij; i_row++)
        {
            for (size_t i_col = 0; i_col != as_size(naux_mu); i_col++)
            {
                (*data)[i_col * n_ij + i_row] = Cs_in[i_row * naux_mu + i_col];
            }
        }
        const std::initializer_list<std::size_t> shape{as_size(naux_mu), as_size(nbasis_i),
                                                       as_size(nbasis_j)};
        cs_data.data_libri[I][{J, Ra}] = RI::Tensor<double>(shape, data);
        cs_data.use_libri = true;
    }
    else
    {
        librpa_int::Vector3_Order<int> box(R[0], R[1], R[2]);
        std::shared_ptr<matrix> cs_ptr = std::make_shared<matrix>();
        cs_ptr->create(nbasis_i * nbasis_j, naux_mu);
        memcpy((*cs_ptr).c, Cs_in, sizeof(double) * cs_size);
        cs_data.data_IJR[I][J][box] = cs_ptr;
        cs_data.use_libri = false;
    }
}

static void _set_aux_coulomb_k_atom_pair(const librpa_int::Vector3_Order<double> &qvec,
                                         librpa_int::atom_t I, librpa_int::atom_t J, size_t naux_mu, size_t naux_nu,
                                         const double* Vq_real_in, const double* Vq_imag_in,
                                         librpa_int::atpair_k_cplx_mat_t &coulomb_mat,
                                         double vq_threshold)
{
    using librpa_int::ComplexMatrix;

    std::shared_ptr<ComplexMatrix> vq_ptr = std::make_shared<ComplexMatrix>();
    vq_ptr->create(naux_mu, naux_nu);

    // Copy real data
    for (size_t i_mu = 0; i_mu < naux_mu; i_mu++)
    {
        for (size_t i_nu = 0; i_nu < naux_nu; i_nu++)
        {
            const auto id = i_nu + i_mu * naux_nu;
            (*vq_ptr)(i_mu, i_nu) = librpa_int::cplxdb(Vq_real_in[id], Vq_imag_in[id]);
        }
    }

    if ((*vq_ptr).real().absmax() >= vq_threshold)
    {
        coulomb_mat[I][J][qvec] = vq_ptr;
    }
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_bare_coulomb_k_atom_pair,
                     int ik, int I, int J, int naux_mu, int naux_nu,
                     const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold)
{
    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &qvec = pds->pbc.klist[ik];
    _set_aux_coulomb_k_atom_pair(qvec, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, pds->vq, vq_threshold);
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_cut_coulomb_k_atom_pair,
                     int ik, int I, int J, int naux_mu, int naux_nu,
                     const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold)
{
    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &qvec = pds->pbc.klist[ik];
    _set_aux_coulomb_k_atom_pair(qvec, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, pds->vq_cut, vq_threshold);
}

static void _set_aux_coulomb_k_2D_block(
    const librpa_int::Vector3_Order<double>& qvec, int mu_begin, int mu_end,
    int nu_begin, int nu_end, const double* Vq_real_in, const double* Vq_imag_in,
    std::map<librpa_int::Vector3_Order<double>, librpa_int::Matz>& vq_block)
{
    using namespace librpa_int;
    using std::shared_ptr;
    using std::make_shared;

    int n_mu_loc = mu_end - mu_begin;
    int n_nu_loc = nu_end - nu_begin;
    if (n_mu_loc < 1 || n_nu_loc < 1) return;

    Matz mat(n_mu_loc, n_nu_loc, MAJOR::ROW);

    size_t ii = 0;
    for (int i_mu = 0; i_mu < n_mu_loc; i_mu++)
    {
        for (int i_nu = 0; i_nu < n_nu_loc; i_nu++)
        {
            mat(i_mu, i_nu) = complex<double>(Vq_real_in[ii], Vq_imag_in[ii]);
            ii += 1;
        }
    }
    mat.swap_to_col_major();
    vq_block[qvec] = std::move(mat);
}

static void _parse_vq_dims(int &lbrow, int &ubrow, int &lbcol, int &ubcol,
                          const int lbrow_in, const int ubrow_in, const int lbcol_in, const int ubcol_in)
{
    if (lbrow < 0 || ubrow < 0 || lbcol < 0 || ubcol < 0)
    {
        lbrow = lbrow_in;
        ubrow = ubrow_in;
        lbcol = lbcol_in;
        ubcol = ubcol_in;
        return;
    }

    if (lbrow != lbrow_in || ubrow != ubrow_in || lbcol != lbcol_in || ubcol != ubcol_in)
    {
        throw LIBRPA_RUNTIME_ERROR("{lb,ub}{row,col} input is inconsistent from previous value");
    }
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_bare_coulomb_k_2d_block,
                     int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                     const double* Vq_real_in, const double* Vq_imag_in)
{
    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &qvec = pds->pbc.klist[ik];
    _parse_vq_dims(pds->vq_lbrow, pds->vq_ubrow, pds->vq_lbcol, pds->vq_ubcol,
                   mu_begin, mu_end, nu_begin, nu_end);
    _set_aux_coulomb_k_2D_block(qvec, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in,
                                Vq_imag_in, pds->vq_block_loc);
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_cut_coulomb_k_2d_block,
                     int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                     const double* Vq_real_in, const double* Vq_imag_in)
{
    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &qvec = pds->pbc.klist[ik];
    _parse_vq_dims(pds->vq_lbrow, pds->vq_ubrow, pds->vq_lbcol, pds->vq_ubcol,
                   mu_begin, mu_end, nu_begin, nu_end);
    _set_aux_coulomb_k_2D_block(qvec, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in,
                                Vq_imag_in, pds->vq_cut_block_loc);
}

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_dielect_func_imagfreq, int nfreq, const double *omegas_imag, const double *dielect_func)
{
    using namespace librpa_int;
    auto pds = librpa_int::api::get_dataset_instance(h);
    pds->omegas_imagfreq = std::vector<double>(omegas_imag, omegas_imag + nfreq);
    pds->epsmacs_imagfreq = std::vector<double>(dielect_func, dielect_func + nfreq);
}
