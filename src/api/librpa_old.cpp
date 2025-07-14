/*
 * @file librpa.cpp
 * @author Rong Shi (srlive@mail.ustc.edu.cn)
 * @brief
 * @version 0.1
 * @date 2023-08-18
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "librpa.h"

#include <stdlib.h>
<<<<<<< HEAD:src/api/librpa_old.cpp

#include <cstring>

#include "../core/atom.h"
#include "../core/meanfield.h"
#include "../core/params.h"
#include "../core/pbc.h"
#include "../core/ri.h"
#include "../math/matrix_m.h"
#include "../mpi/envs_blacs.h"
#include "../mpi/global_mpi.h"
#include "../utils/constants.h"
#include "../io/global_io.h"
#include "librpa_main.h"
=======

#include <cstring>

#include "atoms.h"
#include "constants.h"
#include "envs_blacs.h"
#include "envs_io.h"
#include "envs_mpi.h"
#include "librpa_main.h"
#include "matrix_m.h"
#include "meanfield.h"
#include "parallel_mpi.h"
#include "params.h"
#include "pbc.h"
#include "ri.h"
#include "utils_io.h"
>>>>>>> b1f933c (construct head):src/librpa.cpp

#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>

#include <initializer_list>
#else
#include "../utils/libri_stub.h"
#endif

#include "app_exx.h"
#include "app_rpa.h"

<<<<<<< HEAD:src/api/librpa_old.cpp

void initialize_librpa_environment(
        int redirect_stdout, const char *output_filename)
=======
void initialize_librpa_environment(MPI_Comm comm_global_in, int is_fortran_comm,
                                   int redirect_stdout, const char* output_filename)
>>>>>>> b1f933c (construct head):src/librpa.cpp
{
    librpa_int::global::init_global_mpi();
    librpa_int::envs::initialize_blacs(librpa_int::global::mpi_comm_global);
    librpa_int::global::init_global_io(redirect_stdout, output_filename);

    // cout << mpi_comm_global_h.str() << endl;
}

void finalize_librpa_environment()
{
    librpa_int::global::finalize_global_io();
    librpa_int::envs::finalize_blacs();
    librpa_int::global::finalize_global_mpi();
}

void set_dimension(int nspins, int nkpts, int nstates, int nbasis, int natoms)
{
    using namespace librpa_int;
    using librpa_int::global::mpi_comm_global_h;

    if (mpi_comm_global_h.is_root())
    {
        cout << "In LibRPA nspin: " << nspins << "  nkpt: " << nkpts << endl;
<<<<<<< HEAD:src/api/librpa_old.cpp
        global::lib_printf("In LibRPA nspin: %d,  nkpts: %d nstate: %d, nbasis: %d, natom: %d\n",
                nspins, nkpts, nstates, nbasis, natoms);
=======
        LIBRPA::utils::lib_printf(
            "In LibRPA nspin: %d,  nkpts: %d nstate: %d, nbasis: %d, natom: %d\n", nspins, nkpts,
            nstates, nbasis, natoms);
>>>>>>> b1f933c (construct head):src/librpa.cpp
    }
    meanfield.set(nspins, nkpts, nstates, nbasis);
    natom = natoms;
}

void set_ao_basis_wfc(int ispin, int ik, double* wfc_real, double* wfc_imag, MeanField& mf)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

    meanfield.get_efermi() = efermi;
    auto& eskb = meanfield.get_eigenvals();
    auto& swg = meanfield.get_weight();
=======
    // LIBRPA::utils::lib_printf("is: %d, ik: %d\n",is,ik);
    // int length_ib_iw=meanfield.get_n_bands()*meanfield.get_n_aos();
    // vector<double> vec_wfc_real(wfc_real,wfc_real+length_ib_iw);
    // vector<double> vec_wfc_imag(wfc_imag,wfc_imag+length_ib_iw);
    if (Params::use_soc)
    {
        auto& wfc1 = mf.get_eigenvectors().at(ispin).at(0).at(ik);
        auto& wfc2 = mf.get_eigenvectors().at(ispin).at(1).at(ik);
        const auto nbands = mf.get_n_bands();
        const auto naos = mf.get_n_aos();
        const auto nsoc = mf.get_n_soc();
        const auto nbasis = naos * nsoc;
        const auto ntot = naos * nbands;
        const auto nbs = nbands * nsoc;

        /*for (int isoc = 0; isoc != nsoc; isoc++)
        {
            for (int i = 0; i != naos; i++)
            {
                for (int j = 0; j != nbands; j++)
                {
                    if (isoc == 0)
                        wfc1(j, i) = complex<double>(wfc_real[isoc * ntot + i * nbands + j],
                                                     wfc_imag[isoc * ntot + i * nbands + j]);
                    else if (isoc == 1)
                        wfc2(j, i) = complex<double>(wfc_real[isoc * ntot + i * nbands + j],
                                                     wfc_imag[isoc * ntot + i * nbands + j]);
                }
            }
        }*/

        for (int i = 0; i != naos; i++)
        {
            for (int isoc = 0; isoc != nsoc; isoc++)
            {
                for (int j = 0; j != nbands; j++)
                {
                    if (isoc == 0)
                        wfc1(j, i) = complex<double>(wfc_real[i * nbs + isoc * nbands + j],
                                                     wfc_imag[i * nbs + isoc * nbands + j]);
                    else if (isoc == 1)
                        wfc2(j, i) = complex<double>(wfc_real[i * nbs + isoc * nbands + j],
                                                     wfc_imag[i * nbs + isoc * nbands + j]);
                }
            }
        }
    }
    else
    {
        auto& wfc = mf.get_eigenvectors().at(ispin).at(0).at(ik);
        const auto n = mf.get_n_bands() * mf.get_n_aos();
        for (int i = 0; i != n; i++)
        {
            // LIBRPA::utils::lib_printf("In ao wfc: %f, %f\n",wfc_real[i],wfc_imag[i]);
            wfc.c[i] = complex<double>(wfc_real[i], wfc_imag[i]);
        }
    }

    // print_complex_matrix("wfc_isk", wfc.at(is).at(ik));
}

void set_wg_ekb_efermi(int nspins, int nkpts, int nstates, double* wg, double* ekb, double efermi,
                       MeanField& mf)
{
    mf.get_efermi() = efermi;
    auto& eskb = mf.get_eigenvals();
    auto& swg = mf.get_weight();
>>>>>>> b1f933c (construct head):src/librpa.cpp
    int length_kb = nkpts * nstates;
    for (int is = 0; is != nspins; is++)
    {
        memcpy(eskb[is].c, ekb + length_kb * is, length_kb * sizeof(double));
        memcpy(swg[is].c, wg + length_kb * is, length_kb * sizeof(double));
        // wg[is](k_index, i) = stod(ws) / n_kpoints; // different with abacus!
        swg[is] *= (1.0 / nkpts);
    }
<<<<<<< HEAD:src/api/librpa_old.cpp
    // for(int is=0;is!=nspins;is++)
    // {
    //     print_matrix(" eskb",eskb[is]);
    //     print_matrix(" swg",swg[is]);
    // }
}

void set_ao_basis_wfc(int ispin, int ik, double* wfc_real, double* wfc_imag)
{
    using namespace librpa_int;
    // librpa_int::global::lib_printf("is: %d, ik: %d\n",is,ik);
    // int length_ib_iw=meanfield.get_n_bands()*meanfield.get_n_aos();
    // vector<double> vec_wfc_real(wfc_real,wfc_real+length_ib_iw);
    // vector<double> vec_wfc_imag(wfc_imag,wfc_imag+length_ib_iw);
    auto& wfc = meanfield.get_eigenvectors().at(ispin).at(ik);
    const auto n = meanfield.get_n_bands() * meanfield.get_n_aos();
    for (int i = 0; i != n; i++)
    {
        // librpa_int::global::lib_printf("In ao wfc: %f, %f\n",wfc_real[i],wfc_imag[i]);
        wfc.c[i] = complex<double>(wfc_real[i], wfc_imag[i]);
    }
    // print_complex_matrix("wfc_isk", wfc.at(is).at(ik));
=======
>>>>>>> b1f933c (construct head):src/librpa.cpp
}

void set_latvec_and_G(double lat_mat[9], double G_mat[9])
{
    using namespace librpa_int;

    latvec.e11 = lat_mat[0];
    latvec.e12 = lat_mat[1];
    latvec.e13 = lat_mat[2];

    latvec.e21 = lat_mat[3];
    latvec.e22 = lat_mat[4];
    latvec.e23 = lat_mat[5];

    latvec.e31 = lat_mat[6];
    latvec.e32 = lat_mat[7];
    latvec.e33 = lat_mat[8];

    // latvec /= ANG2BOHR;
    lat_array[0] = {latvec.e11, latvec.e12, latvec.e13};
    lat_array[1] = {latvec.e21, latvec.e22, latvec.e23};
    lat_array[2] = {latvec.e31, latvec.e32, latvec.e33};

    G.e11 = G_mat[0];
    G.e12 = G_mat[1];
    G.e13 = G_mat[2];

    G.e21 = G_mat[3];
    G.e22 = G_mat[4];
    G.e23 = G_mat[5];

    G.e31 = G_mat[6];
    G.e32 = G_mat[7];
    G.e33 = G_mat[8];

    G /= TWO_PI;
    // G *= ANG2BOHR;

    // if (mpi_comm_global_h.is_root())
    // {
<<<<<<< HEAD:src/api/librpa_old.cpp
    //     librpa_int::global::lib_printf(" LibRPA_lat (Bohr) : %f, %f, %f\n",latvec.e11,latvec.e12,latvec.e13);
=======
    //     LIBRPA::utils::lib_printf(" LibRPA_lat (Bohr) : %f, %f,
    //     %f\n",latvec.e11,latvec.e12,latvec.e13);
>>>>>>> b1f933c (construct head):src/librpa.cpp
    // }
    // latvec.print();
    // G.print();
}

void set_kgrids_kvec_tot(int nk1, int nk2, int nk3, double* kvecs)
{
    using namespace librpa_int;

    kv_nmp[0] = nk1;
    kv_nmp[1] = nk2;
    kv_nmp[2] = nk3;

    kvec_c = new Vector3<double>[nk1 * nk2 * nk3];

    for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
    {
        double kx = kvecs[ik * 3];
        double ky = kvecs[ik * 3 + 1];
        double kz = kvecs[ik * 3 + 2];

        kvec_c[ik] = {kx, ky, kz};
        // kvec_c[ik] *= (ANG2BOHR / TWO_PI);
        kvec_c[ik] /= TWO_PI;
        Vector3_Order<double> kvec_tmp(kvec_c[ik]);
        klist.push_back(kvec_tmp);
        kfrac_list.push_back(latvec * kvec_tmp);
        // librpa_int::global::lib_printf("ik: %d, (%f, %f, %f), (%f, %f, %f)\n",
        //         ik, kvec_c[ik].x, kvec_c[ik].y, kvec_c[ik].z,
        //         kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
    }
}

void set_ibz2bz_index_and_weight(const int nk_irk, const int* ibz2bz_index, const double* wk_irk)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

   // librpa_int::global::lib_printf(" nks_irk: %d\n",nk_irk);
   for (int ik_ibz = 0; ik_ibz != nk_irk; ik_ibz++)
   {
       Vector3_Order<double> kvec_ibz = klist[ibz2bz_index[ik_ibz]];
       klist_ibz.push_back(kvec_ibz);
       irk_weight.insert(pair<Vector3_Order<double>, double>(kvec_ibz, wk_irk[ik_ibz]));
       // librpa_int::global::lib_printf("ibz2bz:  %d   kvec_ibz:( %f, %f,
       // %f)\n",ibz2bz_index[ik_ibz],kvec_ibz.x,kvec_ibz.y,kvec_ibz.z);
       // librpa_int::global::lib_printf("irk_weight: %f\n",irk_weight[kvec_ibz]);
       // for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
       // {
       //     if(klist[ik]==kvec_ibz)
       //         map_irk_ks[kvec_ibz].push_back(klist[ik]);
       // }
   }
}

<<<<<<< HEAD:src/api/librpa_old.cpp
void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in, int insert_index_only)
=======
    // LIBRPA::utils::lib_printf(" nks_irk: %d\n",nk_irk);
    for (int ik_ibz = 0; ik_ibz != nk_irk; ik_ibz++)
    {
        Vector3_Order<double> kvec_ibz = klist[ibz2bz_index[ik_ibz]];
        klist_ibz.push_back(kvec_ibz);
        irk_weight.insert(pair<Vector3_Order<double>, double>(kvec_ibz, wk_irk[ik_ibz]));
        // LIBRPA::utils::lib_printf("ibz2bz:  %d   kvec_ibz:( %f, %f,
        // %f)\n",ibz2bz_index[ik_ibz],kvec_ibz.x,kvec_ibz.y,kvec_ibz.z);
        // LIBRPA::utils::lib_printf("irk_weight: %f\n",irk_weight[kvec_ibz]);
        // for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        // {
        //     if(klist[ik]==kvec_ibz)
        //         map_irk_ks[kvec_ibz].push_back(klist[ik]);
        // }
    }
}

void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in,
<<<<<<< HEAD:src/api/librpa_old.cpp
                      int insert_index_only)
>>>>>>> b1f933c (construct head):src/librpa.cpp
{
    using namespace librpa_int;

    atom_nw.insert(pair<atom_t, size_t>(I, nbasis_i));
    atom_mu.insert(pair<atom_t, size_t>(I, naux_mu));
=======
void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in,
                      int insert_index_only, const std::string keyword)
{
=======
                      int insert_index_only, const std::string keyword)
{
>>>>>>> b005208 (enable shrink abfs):src/librpa.cpp
    if (keyword == "Cs_data")
    {
        atom_nw.insert(pair<atom_t, size_t>(I, nbasis_i));
        atom_mu.insert(pair<atom_t, size_t>(I, naux_mu));
    }
<<<<<<< HEAD:src/api/librpa_old.cpp
<<<<<<< HEAD:src/api/librpa_old.cpp
>>>>>>> 7242b2c (enable shrink abfs):src/librpa.cpp
=======
>>>>>>> b005208 (enable shrink abfs):src/librpa.cpp
=======
    else if (keyword == "Cs_shrinked_data")
    {
        atom_mu_s.insert(pair<atom_t, size_t>(I, naux_mu));
    }
>>>>>>> 5b31a70 (fix: read sinvS for multi-mpi):src/librpa.cpp

    if (insert_index_only)
    {
        return;
    }

    // cout<< ia1<<ia2<<box<<endl;
    // librpa_int::global::lib_printf("Cs_in size: %zu",sizeof(Cs_in) / sizeof(Cs_in[0]));
    const int cs_size = nbasis_i * nbasis_j * naux_mu;
    const int n_ij = nbasis_i * nbasis_j;
    //(*cs_ptr).c=Cs_in;

    /*
     * librpa_int::parallel_routing may not be properly set when parsing to set_ao_basis_aux.
     * Therefore using Params::parallel_routing string to check, because
     * Params are required to set up after the environment initialization and before data
     * transfer.
     * */
    // Cs_data.use_libri = librpa_int::parallel_routing == librpa_int::ParallelRouting::LIBRI;
    Cs_data.use_libri = Params::parallel_routing == "libri";
    Cs_shrinked_data.use_libri = Params::parallel_routing == "libri";

    if (Cs_data.use_libri)
    {
        const std::array<int, 3> Ra{R[0], R[1], R[2]};
        // RI tensor uses ABF as slowest index, so we need transpose the input data.
        auto data = make_shared<std::valarray<double>>(cs_size);
        // Unless n_cols >> n_rows, cache-friendly reading (by row in C/C++) is more efficient
        for (int i_row = 0; i_row < n_ij; i_row++)
        {
            for (int i_col = 0; i_col != naux_mu; i_col++)
            {
                (*data)[i_col * n_ij + i_row] = Cs_in[i_row * naux_mu + i_col];
            }
        }
        const std::initializer_list<std::size_t> shape{static_cast<std::size_t>(naux_mu),
                                                       static_cast<std::size_t>(nbasis_i),
                                                       static_cast<std::size_t>(nbasis_j)};
        if (keyword == "Cs_data")
        {
            Cs_data.data_libri[I][{J, Ra}] = RI::Tensor<double>(shape, data);
        }
        else if (keyword == "Cs_shrinked_data")
        {
            Cs_shrinked_data.data_libri[I][{J, Ra}] = RI::Tensor<double>(shape, data);
        }
    }
    else
    {
        Vector3_Order<int> box(R[0], R[1], R[2]);
        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(nbasis_i * nbasis_j, naux_mu);
        memcpy((*cs_ptr).c, Cs_in, sizeof(double) * cs_size);
<<<<<<< HEAD:src/api/librpa_old.cpp
<<<<<<< HEAD:src/api/librpa_old.cpp
        Cs_data.data_IJR[I][J][box] = cs_ptr;
        // librpa_int::global::lib_printf("Cs out:\n");
=======
=======
>>>>>>> b005208 (enable shrink abfs):src/librpa.cpp
        if (keyword == "Cs_data")
        {
            Cs_data.data_IJR[I][J][box] = cs_ptr;
        }
        else if (keyword == "Cs_shrinked_data")
        {
            Cs_shrinked_data.data_IJR[I][J][box] = cs_ptr;
        }
        // LIBRPA::utils::lib_printf("Cs out:\n");
>>>>>>> 7242b2c (enable shrink abfs):src/librpa.cpp
        // for(int i=0;i!=cs_size;i++)
        //     librpa_int::global::lib_printf("   %f",(*Cs[I][J][box]).c[i]);
    }
}

<<<<<<< HEAD:src/api/librpa_old.cpp
static void _set_aux_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu, double* Vq_real_in, double* Vq_imag_in,
                                         librpa_int::atpair_k_cplx_mat_t &coulomb_mat)
=======
static void _set_aux_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu,
                                         double* Vq_real_in, double* Vq_imag_in,
                                         atpair_k_cplx_mat_t& coulomb_mat)
>>>>>>> b1f933c (construct head):src/librpa.cpp
{
    using namespace librpa_int;

    // librpa_int::global::lib_printf("I,J,mu,nu: %d  %d  %d  %d\n",I,J, naux_mu,naux_nu);
    // librpa_int::global::lib_printf("atom_mu nu: %d %d\n",atom_mu[I],atom_mu[J]);
    // librpa_int::global::lib_printf("Vq threshold : %f\n",Params::vq_threshold);
    Vector3_Order<double> qvec(klist[ik]);
    // librpa_int::global::lib_printf("qvec: %f,%f,%f\n",qvec.x,qvec.y,qvec.z);
    shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
    vq_ptr->create(naux_mu, naux_nu);
    // librpa_int::global::lib_printf("vq_ptr_size: %d\n",(*vq_ptr).size);
    // librpa_int::global::lib_printf(" vq_in 0 0 :  %f\n",Vq_real_in[0]);
    for (int i_mu = 0; i_mu != naux_mu; i_mu++)
    {
        for (int i_nu = 0; i_nu != naux_nu; i_nu++)
        {
            // librpa_int::global::lib_printf("mu,nu:  %d %d, vq_real: %f,  vq_imag:
            // %f\n",i_mu,i_nu,Vq_real_in[i_nu+i_mu*naux_nu],Vq_imag_in[i_nu+i_mu*naux_nu]);

            (*vq_ptr)(i_mu, i_nu) = complex<double>(Vq_real_in[i_nu + i_mu * naux_nu],
                                                    Vq_imag_in[i_nu + i_mu * naux_nu]);
        }
    }
    if ((*vq_ptr).real().absmax() >= Params::vq_threshold)
    {
        coulomb_mat[I][J][qvec] = vq_ptr;
    }
    // print_complex_matrix("Vq",(*Vq[I][J][qvec]));
}

void set_aux_bare_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu,
                                      double* Vq_real_in, double* Vq_imag_in)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

    _set_aux_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu,
            Vq_real_in, Vq_imag_in, Vq);
=======
    _set_aux_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, Vq);
>>>>>>> b1f933c (construct head):src/librpa.cpp
}

void set_aux_cut_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu,
                                     double* Vq_real_in, double* Vq_imag_in)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

    _set_aux_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu,
            Vq_real_in, Vq_imag_in, Vq_cut);
}

static void _set_aux_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end,
        double* Vq_real_in, double* Vq_imag_in, std::map<librpa_int::Vector3_Order<double>, librpa_int::ComplexMatrix> &vq_block)
=======
    _set_aux_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, Vq_cut);
}

static void _set_aux_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end,
                                        int nu_begin, int nu_end, double* Vq_real_in,
                                        double* Vq_imag_in,
                                        map<Vector3_Order<double>, ComplexMatrix>& vq_block)
>>>>>>> b1f933c (construct head):src/librpa.cpp
{
    using namespace librpa_int;

    int brow = mu_begin - 1;
    int erow = mu_end - 1;
    int bcol = nu_begin - 1;
    int ecol = nu_end - 1;

    Vector3_Order<double> qvec(klist[ik]);
    // librpa_int::global::lib_printf("qvec: %f,%f,%f\n",qvec.x,qvec.y,qvec.z);
    shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
    vq_ptr->create(max_naux, max_naux);

    if (!Vq_block_loc.count(qvec))
    {
        Vq_block_loc[qvec].create(max_naux, max_naux);
    }
    int ii = 0;
    for (int i_mu = brow; i_mu <= erow; i_mu++)
    {
        for (int i_nu = bcol; i_nu <= ecol; i_nu++)
        {
            vq_block[qvec](i_mu, i_nu) = complex<double>(Vq_real_in[ii], Vq_imag_in[ii]);
            ii += 1;
        }
    }
}

void set_aux_bare_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin,
                                     int nu_end, double* Vq_real_in, double* Vq_imag_in)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

    _set_aux_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in,
            Vq_block_loc);
=======
    _set_aux_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in,
                                Vq_imag_in, Vq_block_loc);
>>>>>>> b1f933c (construct head):src/librpa.cpp
}

void set_aux_cut_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin,
                                    int nu_end, double* Vq_real_in, double* Vq_imag_in)
{
<<<<<<< HEAD:src/api/librpa_old.cpp
    using namespace librpa_int;

    _set_aux_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in,
            Vq_cut_block_loc);
=======
    _set_aux_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in,
                                Vq_imag_in, Vq_cut_block_loc);
>>>>>>> b1f933c (construct head):src/librpa.cpp
}

void set_librpa_params(LibRPAParams* params_c)
{
    using namespace librpa_int;

    Params::output_file = params_c->output_file;
    Params::output_dir = params_c->output_dir;
    Params::tfgrids_type = params_c->tfgrids_type;
    Params::parallel_routing = params_c->parallel_routing;
    Params::DFT_software = params_c->DFT_software;

    Params::nfreq = params_c->nfreq;

    Params::debug = params_c->debug;
    Params::use_scalapack_ecrpa = params_c->use_scalapack_ecrpa;

    Params::gf_R_threshold = params_c->gf_R_threshold;
    Params::cs_threshold = params_c->cs_threshold;
    Params::vq_threshold = params_c->vq_threshold;
    Params::libri_chi0_threshold_C = params_c->libri_chi0_threshold_C;
    Params::libri_chi0_threshold_G = params_c->libri_chi0_threshold_G;
    Params::libri_exx_threshold_C = params_c->libri_exx_threshold_C;
    Params::libri_exx_threshold_D = params_c->libri_exx_threshold_D;
    Params::libri_exx_threshold_V = params_c->libri_exx_threshold_V;
}

void get_default_librpa_params(LibRPAParams* params_c)
{
    // All member of LibRPAParams must be set.
    strcpy(params_c->output_file, "stdout");
    strcpy(params_c->output_dir, "librpa.d");
    strcpy(params_c->parallel_routing, "auto");
    strcpy(params_c->tfgrids_type, "minimax");
    strcpy(params_c->DFT_software, "auto");

    params_c->nfreq = 6;

    params_c->debug = 0;
    params_c->use_scalapack_ecrpa = 0;

    params_c->gf_R_threshold = 0.0e0;
    params_c->cs_threshold = 0.0e0;
    params_c->vq_threshold = 0.0e0;
    params_c->libri_chi0_threshold_C = 0.0e0;
    params_c->libri_chi0_threshold_G = 0.0e0;
    params_c->libri_exx_threshold_C = 0.0e0;
    params_c->libri_exx_threshold_D = 0.0e0;
    params_c->libri_exx_threshold_V = 0.0e0;
    params_c->libri_gw_threshold_C = 0.0e0;
    params_c->libri_gw_threshold_G = 0.0e0;
    params_c->libri_gw_threshold_W = 0.0e0;
}

<<<<<<< HEAD:src/api/librpa_old.cpp

// void run_librpa_main()
// {
//     librpa_int::global::lib_printf("Begin run LibRPA\n");
//     // std::ofstream outputFile("LibRPA_cout.txt");
//     // std::streambuf* originalCoutBuffer = std::cout.rdbuf();
//    // std::cout.rdbuf(outputFile.rdbuf());
//     librpa_main();
//
//     //std::cout.rdbuf(originalCoutBuffer);
//     //outputFile.close();
// }
=======
void run_librpa_main()
{
    LIBRPA::utils::lib_printf("Begin run LibRPA\n");
    // std::ofstream outputFile("LibRPA_cout.txt");
    // std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    // std::cout.rdbuf(outputFile.rdbuf());
    librpa_main();

    // std::cout.rdbuf(originalCoutBuffer);
    // outputFile.close();
}
>>>>>>> b1f933c (construct head):src/librpa.cpp

void get_frequency_grids(int ngrid, double* freqeuncy_grids) {}

void get_rpa_correlation_energy(double* rpa_corr, double* rpa_corr_irk_contrib,
                                std::map<Vector3_Order<double>, ComplexMatrix>& sinvS,
                                const std::string& input_dir, const bool use_shrink_abfs)
{
    using namespace librpa_int;

<<<<<<< HEAD:src/api/librpa_old.cpp
    std::complex<double> rpa_corr_;
    std::vector<std::complex<double>> rpa_corr_irk_contrib_(n_irk_points);

<<<<<<< HEAD:src/api/librpa_old.cpp
    librpa_int::app::get_rpa_correlation_energy_(rpa_corr_, rpa_corr_irk_contrib_);
=======
    // std::cout.rdbuf(originalCoutBuffer);
    // outputFile.close();
}

void get_frequency_grids(int ngrid, double* freqeuncy_grids) {}

void get_rpa_correlation_energy(double* rpa_corr, double* rpa_corr_irk_contrib,
                                std::map<Vector3_Order<double>, ComplexMatrix>& sinvS,
                                const std::string& input_dir, const bool use_shrink_abfs)
{
    std::complex<double> rpa_corr_;
    std::vector<std::complex<double>> rpa_corr_irk_contrib_(n_irk_points);

    LIBRPA::app::get_rpa_correlation_energy_(rpa_corr_, rpa_corr_irk_contrib_, sinvS, input_dir,
                                             use_shrink_abfs);
>>>>>>> 7242b2c (enable shrink abfs):src/librpa.cpp
=======
    LIBRPA::app::get_rpa_correlation_energy_(rpa_corr_, rpa_corr_irk_contrib_, sinvS, input_dir,
                                             use_shrink_abfs);
>>>>>>> b005208 (enable shrink abfs):src/librpa.cpp

    auto dp = reinterpret_cast<double*>(rpa_corr_irk_contrib_.data());

    memcpy(rpa_corr, &rpa_corr_, 2 * sizeof(double));
    memcpy(rpa_corr_irk_contrib, dp, 2 * n_irk_points * sizeof(double));
}

void compute_exx_orbital_energy(int i_state_low, int i_state_high, int n_kpoints_task,
                                const int* i_kpoints_task, double* exx)
{
    const auto exx_vec = librpa_int::app::compute_exx_orbital_energy_(i_state_low, i_state_high,
                                                                  n_kpoints_task, i_kpoints_task);
    memcpy(exx, exx_vec.data(), sizeof(double) * exx_vec.size());
}
