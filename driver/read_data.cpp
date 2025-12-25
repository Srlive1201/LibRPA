#include "read_data.h"

#include <dirent.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "driver.h"
#include "librpa.hpp"
#include "../src/mpi/global_mpi.h"
// #include "../src/utils/constants.h"
#include "../src/api/instance_manager.h"
#include "../src/io/global_io.h"
#include "../src/io/stl_io_helper.h"
#include "../src/utils/profiler.h"
#include "../src/utils/utils_mem.h"

using std::ifstream;
using std::string;

void read_scf_occ_eigenvalues(const string &file_path)
{
    using std::to_string;
    using driver::n_spins;
    using driver::n_kpoints;
    using driver::n_states;
    using driver::n_basis_wfc;
    using driver::iks_eigvec_local;
    using librpa_int::global::myid_global;
    using librpa_int::global::size_global;

    // cout << "Begin to read aims-band_out" << endl;
    ifstream infile;
    infile.open(file_path);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    string ks, ss, a, ws, es, d;
    double efermi;
    infile >> n_kpoints;
    infile >> n_spins;
    infile >> n_states;
    infile >> n_basis_wfc;
    infile >> efermi;

    iks_eigvec_local.clear();
    if (driver::get_bool(driver::opts.use_kpara_scf_eigvec))
    {
        for (int ik = 0; ik < driver::n_kpoints; ik++)
        {
            if (ik % size_global == myid_global) iks_eigvec_local.emplace_back(ik);
        }
    }
    else
    {
        for (int ik = 0; ik < driver::n_kpoints; ik++)
            iks_eigvec_local.emplace_back(ik);
    }

    driver::h.set_scf_dimension(n_spins, n_kpoints, n_states, n_basis_wfc,
                                0, n_states, 0, n_basis_wfc);
    driver::n_ibz_kpoints = n_kpoints;

    // Load the file data
    auto eskb = new double [n_spins * n_kpoints * n_states];
    auto wskb = new double [n_spins * n_kpoints * n_states];

    const int n_kb = n_kpoints * n_states;

    int iline = 6;

    //cout<<"|eskb: "<<endl;
    for (int ik = 0; ik != n_kpoints; ik++)
    {
        for (int is = 0; is != n_spins; is++)
        {
            infile >> ks >> ss;
            if (!infile.good())
            {
                throw std::logic_error("Error in reading k- and spin- index: line " + to_string(iline) +
                                       ", file: " + file_path);
            }
            iline++;
            //cout<<ik<<is<<endl;
            int k_index = stoi(ks) - 1;
            // int s_index = stoi(ss) - 1;
            for (int i = 0; i != n_states; i++)
            {
                // iband weight energy(Ha) energy(eV)
                infile >> a >> ws >> es >> d;
                if (!infile.good())
                {
                    throw std::logic_error("Error in reading band energy and occupation: line " + to_string(iline) +
                                           ", file: " + file_path);
                }
                iline++;
                wskb[is * n_kb + k_index * n_states + i] = stod(ws); // different with abacus!
                eskb[is * n_kb + k_index * n_states + i] = stod(es);
                //cout<<" i_band: "<<i<<"    eskb: "<<eskb[is](k_index, i)<<endl;
            }
        }
    }
    // for (int is = 0; is != n_spins; is++)
    //     print_matrix("eskb_mat",eskb[is]);

    driver::h.set_wg_ekb_efermi(n_spins, n_kpoints, n_states, wskb, eskb, efermi);

    // free buffer
    delete [] eskb;
    delete [] wskb;
}

int read_vxc(const string &file_path, std::vector<matrix> &vxc)
{
    ifstream infile;
    infile.open(file_path);
    double ha, ev;
    int n_spins, n_kpoints, n_states;
    // int retcode;

    // dimension information
    infile >> n_kpoints;
    infile >> n_spins;
    infile >> n_states;
    if (!infile.good())
    {
        return 1;
    }

    vxc.clear();
    vxc.resize(n_spins);
    for (int is = 0; is != n_spins; is++)
    {
        vxc[is].create(n_kpoints, n_states);
    }

    for (int ik = 0; ik != n_kpoints; ik++)
    {
        for (int is = 0; is != n_spins; is++)
        {
            for (int i = 0; i != n_states; i++)
            {
                infile >> ha >> ev;
                if (!infile.good())
                {
                    return 2;
                }
                vxc[is](ik, i) = ha;
            }
        }
    }
    return 0;
}

static int handle_KS_file(const string &file_path)
{
    using driver::iks_eigvec_local;

    int ret = 0;
    // cout<<file_path<<endl;
    ifstream infile;
    // cout << "Reading eigenvector from file " << file_path << endl;
    infile.open(file_path);
    if (!infile.good())
        return 1;

    string rvalue, ivalue, kstr;

    const auto nspin = driver::n_spins;
    const auto nband = driver::n_states;
    const auto nao = driver::n_basis_wfc;
    const auto n = nband * nao;

    std::vector<double> re(nspin * nband * nao);
    std::vector<double> im(nspin * nband * nao);

    while (infile.peek() != EOF)
    {
        infile >> kstr;
        int ik = stoi(kstr) - 1;
        // cout<<"     ik: "<<ik<<endl;
        if (infile.peek() == EOF)
            break;
        // for aims !!!
        bool skip_this_ik = false;
        if (driver::get_bool(driver::opts.use_kpara_scf_eigvec))
        {
            const auto it = std::find(iks_eigvec_local.cbegin(), iks_eigvec_local.cend(), ik);
            // this k does not belong to this process
            skip_this_ik = (it == iks_eigvec_local.cend());
        }
        for (int iw = 0; iw != nao; iw++)
        {
            for (int ib = 0; ib != nband; ib++)
            {
                for (int is = 0; is != nspin; is++)
                {
                    // cout<<iw<<ib<<is<<ik;
                    infile >> rvalue >> ivalue;
                    if (infile.bad())
                    {
                        ret = 1;
                        break;
                    }
                    if (!skip_this_ik)
                    {
                        // cout<<rvalue<<ivalue<<endl;
                        re[is * n + ib * nao + iw] = stod(rvalue);
                        im[is * n + ib * nao + iw] = stod(ivalue);
                    }
                }
            }
        }
        if (skip_this_ik) continue;
        for (int is = 0; is != nspin; is++)
        {
            driver::h.set_wfc(is, ik, driver::n_states, driver::n_basis_wfc, re.data() + is * n, im.data() + is * n);
        }
        // for abacus
        // for (int ib = 0; ib != NBANDS; ib++)
        //     for (int iw = 0; iw != NLOCAL; iw++)
        //         for (int is = 0; is != NSPIN; is++)
        //         {
        //             // cout<<iw<<ib<<is<<ik;
        //             infile >> rvalue >> ivalue;
        //             // cout<<rvalue<<ivalue<<endl;
        //             wfc_k.at(stoi(ik) - 1)(ib, iw) = complex<double>(stod(rvalue), stod(ivalue));
        //         }
    }
    return ret;
}

int read_eigenvector(const string &dir_path)
{
    // return code
    int ret = 0;
    int files_read = 0;

    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    std::vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        // cout << fm << " find:" << fm.find("KS_eigenvector") << "\n";
        if (fm.find("KS_eigenvector") == 0)
        {
            ret = handle_KS_file(dir_path + fm);
            if (ret != 0)
            {
                break;
            }
            files_read++;
        }
    }
    closedir(dir);
    dir = NULL;

    if (files_read == 0)
    {
        ret = -1;
    }

    //auto tmp_wfc=mf.get_eigenvectors();
    // for(int is=0;is!=mf.get_n_spins();is++)
    //     print_complex_matrix("wfc ",tmp_wfc.at(is).at(0));
    // cout << "Finish read KS_eignvector! " << endl;
    return ret;
}

void read_ri(const string &dir_path, librpa::ParallelRouting &routing)
{
    using driver::n_atoms;
    using driver::n_kpoints;
    using driver::local_atpair;
    using librpa_int::generate_atom_pair_from_nat;
    using librpa_int::decide_auto_routing;
    using librpa_int::dispatch_upper_triangular_tasks;
    using namespace librpa_int::global;

    mpi_comm_global_h.barrier();
    lib_printf_root("Loading RI file from directory: %s\n", dir_path.c_str());

    const auto tot_atpair = generate_atom_pair_from_nat(n_atoms, false);
    const auto tot_atpair_ordered = generate_atom_pair_from_nat(n_atoms, true);

    if (routing == librpa::ParallelRouting::AUTO)
    {
        routing = decide_auto_routing(n_atoms, driver::opts.nfreq * n_kpoints);
    }

    auto pds = librpa_int::api::get_dataset_instance(driver::h.get_c_handler());
    const auto &Cs_data = pds->cs_data;
    const auto &blacs_h = pds->blacs_h;

    local_atpair.clear();

    // HACK: local_atpair should be set in the same mechanism as inside the dataset object,
    //       which is implemented in initialize_ds_atpairs_local in dataset_helper.cpp.
    //       It consists of distributed atom pairs of only upper half, since repsonse function matrix is Hermitian.
    if(routing == librpa::ParallelRouting::ATOMPAIR)
    {
        lib_printf_root("Triangular dispatching of atom pairs\n");
        auto tri_local_atpair = librpa_int::dispatch_upper_triangular_tasks(
            n_atoms, blacs_h.myid, blacs_h.nprows, blacs_h.npcols,
            blacs_h.myprow, blacs_h.mypcol);
        for (const auto &p: tri_local_atpair)
            local_atpair.push_back(p);
        profiler.start("driver_read_Cs");
        read_Cs(dir_path, driver::opts.cs_threshold, local_atpair);
        profiler.stop("driver_read_Cs");
        profiler.start("driver_read_Vq");
        read_Vq_row(dir_path, "coulomb_mat", driver::opts.vq_threshold, local_atpair, false);
        profiler.stop("driver_read_Vq");
    }
    else if(routing == librpa::ParallelRouting::LIBRI)
    {
        lib_printf_root("Evenly distributed Cs and V for LibRI\n");
        profiler.start("driver_read_Cs");
        read_Cs_evenly_distribute(dir_path, driver::opts.cs_threshold, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs);
        profiler.stop("driver_read_Cs");

        lib_printf_coll("| Process %5d: Cs with %14zu non-zero keys from local atpair size %7zu. "
                        "Data memory: %10.2f MB. Wall/CPU time [min]: %12.4f %12.4f\n",
                        mpi_comm_global_h.myid, Cs_data.n_keys(), local_atpair.size(),
                        Cs_data.n_data_bytes() * 8.0e-6,
                        librpa_int::global::profiler.get_wall_time_last("driver_read_Cs") / 60.0,
                        librpa_int::global::profiler.get_cpu_time_last("driver_read_Cs") / 60.0);
        // Vq distributed using the same strategy
        // There should be no duplicate for V

        librpa_int::global::profiler.start("driver_read_Vq");
        auto trangular_loc_atpair = librpa_int::dispatch_upper_triangular_tasks(
            n_atoms, blacs_h.myid, blacs_h.nprows, blacs_h.npcols,
            blacs_h.myprow, blacs_h.mypcol);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        read_Vq_row(dir_path, "coulomb_mat", driver::opts.vq_threshold, local_atpair, false);
        mpi_comm_global_h.barrier();
        librpa_int::global::profiler.stop("driver_read_Vq");
        lib_printf_coll("| Process %5d: coulomb_mat read. Wall/CPU time [min]: %12.4f %12.4f\n",
                        mpi_comm_global_h.myid,
                        librpa_int::global::profiler.get_wall_time_last("driver_read_Vq") / 60.0,
                        librpa_int::global::profiler.get_cpu_time_last("driver_read_Vq") / 60.0);
    }
    else
    {
        lib_printf_root("Complete copy of Cs and V on each process\n");
        local_atpair = generate_atom_pair_from_nat(n_atoms, false);
        profiler.start("driver_read_Cs");
        read_Cs(dir_path, driver::opts.cs_threshold, local_atpair);
        profiler.stop("driver_read_Cs");

        profiler.start("driver_read_Vq");
        read_Vq_full(dir_path, "coulomb_mat", false);
        profiler.stop("driver_read_Vq");
    }
}

//! Check if Cs data file is in ASCII text or unformatted binary format
static bool check_Cs_file_binary(const string &file_path)
{
    // Current strategy:
    //   Assume the file is ASCII, try to read to the first integer, which is the number of atoms
    //   If it succeeds, then the file is ASCII, otherwise it is unformatted.
    //
    // This is the simplest way, and gives less false positives than assuming the file is binary
    // and checking the first integer by reading the first 4 bytes with infile.read.
    bool is_binary = true;
    ifstream infile;
    int natom;
    // infile.open(file_path, std::ios::in | std::ios::binary);
    infile.open(file_path, std::ios::in);
    // infile.read((char *) &natom, sizeof(int));
    infile >> natom;
    if (infile.good())
    {
        is_binary = false;
    }
    // cout << natom << " " << is_binary << endl;
    infile.close();
    return is_binary;
}

//! Check if Coulomb matrix data file is in ASCII text or unformatted binary format
static bool check_coulomb_file_binary(const string &file_path)
{
    bool is_binary = true;
    ifstream infile;
    int nirk;
    infile.open(file_path, std::ios::in);
    infile >> nirk;
    if (infile.good())
    {
        is_binary = false;
    }
    // cout << nirk << " " << is_binary << endl;
    infile.close();
    return is_binary;
}

static size_t handle_Cs_file(const string &file_path, double threshold, const std::vector<atpair_t> &local_atpair)
{
    using namespace std;

    set<size_t> loc_atp_index;
    for(auto &lap:local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    // cout<<"READING Cs from file: "<<file_path<<"  Cs_first_size: "<<loc_atp_index.size()<<endl;
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    size_t cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    int R[3];
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // int natom = stoi(natom_s);
    // int ncell = stoi(ncell_s);

    /* cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl; */
    // for(int loop=0;loop!=natom*natom*ncell;loop++)
    while (infile.peek() != EOF)
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF)
            break;
        // cout << " ia1_s,ia2_s: " << ia1_s << "  " << ia2_s << endl;
        infile >> j_s >> mu_s;
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = stoi(ia1_s) - 1;
        int ia2 = stoi(ia2_s) - 1;
        R[0] = stoi(ic_1);
        R[1] = stoi(ic_2);
        R[2] = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        // cout<< ia1<<ia2<<box<<endl;
        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(n_i * n_j, n_mu);
        // cout<<cs_ptr->nr<<cs_ptr->nc<<endl;

        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> Cs_ele;
                    (*cs_ptr)(i * n_j + j, mu) = stod(Cs_ele);
                    // if (i == j)
                    // {
                    //     (*cs_ptr)(i * n_j + j, mu) = 1.0;
                    // }
                }
        // if(!loc_atp_index.count(ia1))
        //     continue;
        // if (box == Vector3_Order<int>({0, 0, 1}))continue;
        bool keep = loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold;
        if (keep)
            driver::h.set_lri_coeff(driver::opts.parallel_routing, ia1, ia2, n_i, n_j, n_mu, R, cs_ptr->c);
        // cout<<cs_ptr->nr<<cs_ptr->nc<<endl;
        if (!keep)
        {
            cs_discard++;
        }
    }
    infile.close();
    return cs_discard;
}

static size_t handle_Cs_file_binary(const string &file_path, double threshold, const std::vector<atpair_t> &local_atpair)
{
    using namespace std;

    set<size_t> loc_atp_index;
    for(auto &lap:local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    // cout<<"READING Cs from file: "<<file_path<<"  Cs_first_size: "<<loc_atp_index.size()<<endl;
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    size_t cs_discard = 0;
    ifstream infile;
    int dims[8];
    int n_apcell_file;
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));

    int R[3];

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        R[0] = dims[2];
        R[1] = dims[3];
        R[2] = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        // cout<< ia1<<ia2<<box<<endl;

        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(n_i * n_j, n_mu);
        infile.read((char *) cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
        bool keep = loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold;
        // cout << (*cs_ptr).absmax() << "\n";
        if (keep)
        {
            driver::h.set_lri_coeff(driver::opts.parallel_routing, ia1, ia2, n_i, n_j, n_mu, R,
                                    cs_ptr->c);
        }
        else
        {
            cs_discard++;
        }
    }
    return cs_discard;
}

size_t read_Cs(const string &dir_path, double threshold, const std::vector<atpair_t> &local_atpair)
{
    using namespace std;

    size_t cs_discard = 0;
    // cout << "Begin to read Cs" << endl;
    // cout << "cs_threshold:  " << threshold << endl;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    std::vector<string> files;
    bool binary;
    bool binary_checked = false;

    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("Cs_data") == 0)
        {
            const auto fn = dir_path + fm;
            if (!binary_checked)
            {
                binary = check_Cs_file_binary(fn);
                binary_checked = true;
                if (librpa_int::global::myid_global == 0)
                {
                    if (binary)
                    {
                        cout << "Unformatted binary Cs files detected" << endl;
                    }
                    else
                    {
                        cout << "ASCII format Cs files detected" << endl;
                    }
                }
            }
            if (binary)
            {
                cs_discard += handle_Cs_file_binary(fn, threshold, local_atpair);
            }
            else
            {
                cs_discard += handle_Cs_file(fn, threshold, local_atpair);
            }
        }
    }
    closedir(dir);
    dir = NULL;
    // initialize basis set object
    // librpa_int::atomic_basis_wfc.set(atom_nw);
    // librpa_int::atomic_basis_abf.set(atom_mu);
    
    // atom_mu_part_range.resize(atom_mu.size());
    // atom_mu_part_range[0]=0;
    // for(int I=1;I!=atom_mu.size();I++)
    //     atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    // N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // init_N_all_mu(); // FIXME: backward compat

    // for(int i=0;i!=atom_mu_part_range.size();i++)
    //     cout<<" atom_mu_part_range ,i: "<<i<<"    "<<atom_mu_part_range[i]<<endl;

    // cout << "Finish read Cs" << endl;
    return cs_discard;
}

std::vector<size_t> handle_Cs_file_dry(const string &file_path, double threshold)
{
    using namespace std;
    using namespace librpa_int;

    std::vector<size_t> Cs_ids_keep;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // int natom = stoi(natom_s);
    // int ncell = stoi(ncell_s);

    size_t id = 0;
    // int R[3];

    while (infile.peek() != EOF)
    {
        infile >> ia1_s;
        if (infile.peek() == EOF)
            break;
        infile >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s >> j_s >> mu_s;
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);
        // int ia1 = stoi(ia1_s);
        // int ia2 = stoi(ia2_s);
        // R[0] = stoi(ic_1);
        // R[1] = stoi(ic_2);
        // R[2] = stoi(ic_3);

        double maxval = -1.0;
        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> Cs_ele;
                    maxval = std::max(maxval, std::abs(stod(Cs_ele)));
                }
        librpa_int::global::ofs_myid << id << " (" << ic_1 << "," << ic_2 << "," << ic_3 << ") " << maxval << " keep? " << (maxval >= threshold) << endl;
        if (maxval >= threshold)
            Cs_ids_keep.push_back(id);
        id++;
    }
    librpa_int::global::ofs_myid << file_path << ": " << Cs_ids_keep << endl;
    infile.close();
    return Cs_ids_keep;
}

std::vector<size_t> handle_Cs_file_binary_dry(const string &file_path, double threshold)
{
    std::vector<size_t> Cs_ids_keep;
    ifstream infile;
    int dims[8];
    int n_apcell_file;
    // int n_processed = 0;
    // int R[3];
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));

    for (int i_file = 0; i_file < n_apcell_file; i_file++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        // const int ia1 = dims[0] - 1;
        // const int ia2 = dims[1] - 1;
        // R[0] = dims[2];
        // R[1] = dims[3];
        // R[2] = dims[4];
        const int n_i = dims[5];
        const int n_j = dims[6];
        const int n_mu = dims[7];
        // set_ao_basis_aux(ia1, ia2, n_i, n_j, n_mu, R, nullptr, 1);

        matrix mat(n_i * n_j, n_mu);
        infile.read((char *) mat.c, n_i * n_j * n_mu * sizeof(double));
        double maxval = mat.absmax();
        // n_processed++;
        if (maxval >= threshold)
        {
            Cs_ids_keep.push_back(i_file);
#ifdef LIBRPA_DEBUG
            // librpa_int::envs::ofs_myid << i_file << " (" << ic1 << "," << ic2 << "," << ic3 << ") " << maxval << " kept, maxval: " << maxval << endl;
#endif
        }
    }
    // librpa_int::envs::ofs_myid << file_path << ": kept " << Cs_ids_keep.size() << " of " << n_processed << endl;
#ifdef LIBRPA_DEBUG
    // librpa_int::envs::ofs_myid << Cs_ids_keep << endl;
#endif
    infile.close();
    return Cs_ids_keep;
}

static size_t handle_Cs_file_by_ids(const string &file_path, double threshold, const std::vector<size_t> &ids)
{
    using namespace std;
    size_t cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    // int natom, ncell;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // natom = stoi(natom_s);
    // ncell = stoi(ncell_s);
    /* cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl; */
    // for(int loop=0;loop!=natom*natom*ncell;loop++)
    size_t id = 0;
    int R[3];

    while (infile.peek() != EOF)
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF)
            break;
        // cout << " ia1_s,ia2_s: " << ia1_s << "  " << ia2_s << endl;
        infile >> j_s >> mu_s;
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = stoi(ia1_s) - 1;
        int ia2 = stoi(ia2_s) - 1;
        R[0] = stoi(ic_1);
        R[1] = stoi(ic_2);
        R[2] = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        if (std::find(ids.cbegin(), ids.cend(), id) != ids.cend())
        {
            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);

            for (int i = 0; i != n_i; i++)
                for (int j = 0; j != n_j; j++)
                    for (int mu = 0; mu != n_mu; mu++)
                    {
                        infile >> Cs_ele;
                        (*cs_ptr)(i * n_j + j, mu) = stod(Cs_ele);
                    }
            driver::h.set_lri_coeff(driver::opts.parallel_routing, ia1, ia2, n_i, n_j, n_mu, R, cs_ptr->c);
        }
        else
        {
            // set_ao_basis_aux(ia1, ia2, n_i, n_j, n_mu, R, nullptr, 1);

            double maxval = -1.0;
            for (int i = 0; i != n_i; i++)
                for (int j = 0; j != n_j; j++)
                    for (int mu = 0; mu != n_mu; mu++)
                    {
                        infile >> Cs_ele;
                        maxval = std::max(maxval, std::abs(stod(Cs_ele)));
                    }
            if (maxval < threshold) cs_discard++;
        }
        id++;
    }
    infile.close();
    return cs_discard;
}


static size_t handle_Cs_file_binary_by_ids(const string &file_path, double threshold, const std::vector<size_t> &ids)
{
    using namespace std;

    ifstream infile;
    int dims[8];
    int n_apcell_file;
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));
    size_t cs_discard = 0;

    int R[3];

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        R[0] = dims[2];
        R[1] = dims[3];
        R[2] = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        if (std::find(ids.cbegin(), ids.cend(), static_cast<size_t>(i)) != ids.cend())
        {
            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);
            infile.read((char *) cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
            driver::h.set_lri_coeff(driver::opts.parallel_routing, ia1, ia2, n_i, n_j, n_mu, R, cs_ptr->c);
        }
        else
        {
            infile.seekg(n_i * n_j * n_mu * sizeof(double), ios::cur);
            cs_discard++;
        }
    }
    infile.close();
    return cs_discard;
}


size_t read_Cs_evenly_distribute(const string &dir_path, double threshold, int myid, int nprocs)
{
    using namespace std;
    using namespace librpa_int;
    using namespace librpa_int::global;

    size_t cs_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    std::vector<string> files;
    unordered_map<string, std::vector<size_t>> files_Cs_ids;
    unordered_map<string, std::vector<size_t>> files_Cs_ids_this_proc;
    bool binary;
    bool binary_checked = false;

    profiler.start("handle_Cs_file_dry");
    while ((ptr = readdir(dir)) != NULL)
    {
        string fn(ptr->d_name);
        if (fn.find("Cs_data") == 0)
        {
            files.push_back(dir_path + fn);
            if (!binary_checked)
            {
                binary = check_Cs_file_binary(dir_path + fn);
                binary_checked = true;
            }
        }
    }

    const int nfiles = files.size();
    // cout << nfiles << "\n";

    // TODO: the IO can be improved, in two possible ways
    // 1. Each MPI task reads only a subset of files, instead of all files.
    // 2. Parallel reading for each file. This may be more efficient, but would be more difficult to implement
    for (int i_fn = 0; i_fn != nfiles; i_fn++)
    {
        // Let each MPI process read different files at one time
        auto i_fn_myid = (i_fn + myid * nfiles / nprocs) % files.size();
        const auto &fn = files[i_fn_myid];
        ofs_myid << "Reading " << fn << endl;
        std::vector<size_t> ids_keep_this_file;
        if (binary)
        {
            ids_keep_this_file = handle_Cs_file_binary_dry(fn, threshold);
        }
        else
        {
            ids_keep_this_file = handle_Cs_file_dry(fn, threshold);
        }
        files_Cs_ids[fn] = ids_keep_this_file;
    }

    // Filter out the Cs to be actually read in each process
    int id_total = 0;
    for (int i_fn = 0; i_fn < nfiles; i_fn++)
    {
        const auto &fn = files[i_fn];
        const auto &ids_this_file = files_Cs_ids[fn];
        const int n_ids = ids_this_file.size();
        for (int id = 0; id != n_ids; id++)
        {
            if (id_total % nprocs == myid) files_Cs_ids_this_proc[fn].push_back(ids_this_file[id]);
            id_total++;
        }
    }
    profiler.stop("handle_Cs_file_dry");
    closedir(dir);
    dir = NULL;
    if (myid == 0) librpa_int::global::lib_printf("Finished Cs filtering\n");

    profiler.start("handle_Cs_file");
    // cout << files_Cs_ids_this_proc.size() << "\n";
    ofs_myid << "Number of Cs files to process: " << files_Cs_ids_this_proc.size() << "\n";
    for (const auto& fn_ids: files_Cs_ids_this_proc)
    {
        ofs_myid << fn_ids.first << " " << fn_ids.second << endl;
        if (binary)
        {
            cs_discard += handle_Cs_file_binary_by_ids(fn_ids.first, threshold, fn_ids.second);
        }
        else
        {
            cs_discard += handle_Cs_file_by_ids(fn_ids.first, threshold, fn_ids.second);
        }
    }
    librpa_int::global::profiler.stop("handle_Cs_file");

    // initialize basis set object
    // librpa_int::atomic_basis_wfc.set(atom_nw);
    // librpa_int::atomic_basis_abf.set(atom_mu);
    
    // atom_mu_part_range.resize(atom_mu.size());
    // atom_mu_part_range[0]=0;
    // for(int I=1;I!=atom_mu.size();I++)
    //     atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    //
    // N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // cout << "Done\n";
    return cs_discard;
}

void get_natom_ncell_from_first_Cs_file(int &n_atom, int &n_cell, const string &dir_path)
{
    using namespace std;

    // cout<<file_path<<endl;
    ifstream infile;
    bool binary;

    string file_path = "";

    // find Cs file
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    while ((ptr = readdir(dir)) != NULL)
    {
        string fn(ptr->d_name);
        if (fn.find("Cs_data") == 0)
        {
            file_path = dir_path + fn;
            break;
        }
    }
    if (file_path == "")
        throw std::runtime_error("Cs_data file is not found under dir_path: " + dir_path);

    binary = check_Cs_file_binary(file_path);
    if (librpa_int::global::myid_global == 0)
    {
        if (binary)
        {
            cout << "Unformatted binary Cs files detected" << endl;
        }
        else
        {
            cout << "ASCII format Cs files detected" << endl;
        }
    }

    if (binary)
    {
        infile.open(file_path, std::ios::in | std::ios::binary);
        infile.read((char *) &n_atom, sizeof(int));
        infile.read((char *) &n_cell, sizeof(int));
        infile.close();
    }
    else
    {
        string natom_s, ncell_s;
        infile.open(file_path);
        infile >> natom_s >> ncell_s;
        // cout<<"  natom_s:"<<natom_s<<"  ncell_s: "<<ncell_s<<endl;
        n_atom = stoi(natom_s);
        n_cell = stoi(ncell_s);
        infile.close();
    }
}

void read_dielec_func(const string &file_path, std::vector<double> &omegas, std::vector<double> &dielec_func_imagfreq)
{
    std::ifstream ifs;
    double omega, re, im;
    ifs.open(file_path);

    if (!ifs.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    while(ifs >> omega >> re >> im)
    {
        omegas.push_back(omega);
        dielec_func_imagfreq.push_back(re);
    }
    ifs.close();
}

// TODO work-in-progress: 2025-11-30
void read_coulomb(const string &dir_path, const librpa::ParallelRouting routing, bool is_cut)
{
}

static int handle_Vq_full_file(const string &file_path, std::map<int, librpa_int::ComplexMatrix> &Vq_full, bool binary)
{
    // cout << "Begin to read aims vq_real from " << file_path << endl;
    ifstream infile;
    int n_irk_points_local;
    int n_irk_points;

    if (binary)
    {
        infile.open(file_path, std::ios::in | std::ios::binary);
        infile.read((char *) &n_irk_points, sizeof(int));
        infile.read((char *) &n_irk_points_local, sizeof(int));
    }
    else
    {
        infile.open(file_path);
        infile >> n_irk_points;
    }

    if (!infile.good())
        return 1;

    if (binary)
    {
        int nbasbas, brow, erow, bcol, ecol, iq;
        double q_weight;

        for (int i_irk = 0; i_irk < n_irk_points_local; i_irk++)
        {
            infile.read((char *) &nbasbas, sizeof(int));
            infile.read((char *) &brow, sizeof(int));
            infile.read((char *) &erow, sizeof(int));
            infile.read((char *) &bcol, sizeof(int));
            infile.read((char *) &ecol, sizeof(int));
            infile.read((char *) &iq, sizeof(int));
            infile.read((char *) &q_weight, sizeof(double));

            brow--;
            erow--;
            bcol--;
            ecol--;
            iq--;

            if (!Vq_full.count(iq))
            {
                Vq_full[iq].create(nbasbas, nbasbas);
            }

            const int nrow = erow - brow + 1;
            const int ncol = ecol - bcol + 1;
            const size_t n = nrow * ncol;
            std::vector<std::complex<double>> tmp(n);
            infile.read((char *) tmp.data(), 2 * n * sizeof(double));
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    const auto i_mu = i + brow;
                    const auto i_nu = j + bcol;
                    Vq_full[iq](i_mu, i_nu) = tmp[i * ncol + j]; // for abacus
                }
            }
        }
    }
    else
    {
        string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3, vq_r, vq_i, q_num, q_weight;
        while (infile.peek() != EOF)
        {
            infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col;
            if (infile.peek() == EOF)
                break;
            if (!infile.good())
                return 2;
            //cout << "vq range: " << begin_row << " ~ " << end_row << "  ,   " << begin_col << " ~ " << end_col << endl;
            infile >> q_num >> q_weight;
            if (!infile.good())
                return 3;
            int mu = stoi(nbasbas);
            int nu = stoi(nbasbas);
            int brow = stoi(begin_row) - 1;
            int erow = stoi(end_row) - 1;
            int bcol = stoi(begin_col) - 1;
            int ecol = stoi(end_col) - 1;
            int iq = stoi(q_num) - 1;

            //skip empty coulumb_file
            if((erow-brow<=0) || (ecol-bcol<=0) || iq<0)
                return 4;

            if (!Vq_full.count(iq))
            {
                Vq_full[iq].create(mu, nu);
            }
            for (int i_mu = brow; i_mu <= erow; i_mu++)
            {
                for (int i_nu = bcol; i_nu <= ecol; i_nu++)
                {
                    infile >> vq_r >> vq_i;
                    //Vq_full[qvec](i_nu, i_mu) = complex<double>(stod(vq_r), stod(vq_i)); // for FHI-aims
                    Vq_full[iq](i_mu, i_nu) = std::complex<double>(stod(vq_r), stod(vq_i)); // for abacus
                }
            }
        }
    }
    return 0;
}

size_t read_Vq_full(const string &dir_path, const string &vq_fprefix, bool is_cut_coulomb)
{
    using std::cout;
    using std::endl;
    using librpa_int::ComplexMatrix;
    using namespace librpa_int::global;

    auto ds = librpa_int::api::get_dataset_instance(driver::h.get_c_handler());
    const auto &basis_aux = ds->basis_aux;
    const auto atom_mu_part_range = basis_aux.get_part_range();

    size_t vq_save = 0;
    size_t vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    std::vector<string> files;
    std::map<int, ComplexMatrix> Vq_full;

    bool binary;
    bool binary_checked = false;

    profiler.start("handle_Vq_full_file");
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find(vq_fprefix) == 0)
        {
            string file_path = dir_path + fm;
            if (!binary_checked)
            {
                binary = check_coulomb_file_binary(file_path);
                binary_checked = true;
                if (librpa_int::global::myid_global == 0)
                {
                    if (binary)
                    {
                        cout << "Unformatted binary V files detected" << endl;
                    }
                    else
                    {
                        cout << "ASCII format V files detected" << endl;
                    }
                }
            }

            int retcode = handle_Vq_full_file(file_path, Vq_full, binary);
            if (retcode != 0)
            {
                librpa_int::global::lib_printf("Error encountered when reading %s, return code %d", fm.c_str(), retcode);
            }
        }
    }
    profiler.stop("handle_Vq_full_file");

    // cout << "FINISH coulomb files reading!" << endl;
    profiler.start("set_aux_coulomb_k_atom_pair_out");
    for (auto &vf_p : Vq_full)
    {
        int iq = vf_p.first;

        // cout << "Qvec:" << qvec << endl;
        for (size_t I = 0; I != basis_aux.n_atoms; I++)
        {
            for (size_t J = 0; J != basis_aux.n_atoms; J++)
            {
                // Coulomb is Hermitian, only parse upper half
                if (I > J)
                {
                    continue;
                }

                // Vq_full stores the full matrix, parse by I-J block
                // The matrices have to be duplicated ...
                matrix re(basis_aux[I], basis_aux[J]), im(basis_aux[I], basis_aux[J]);

                // vq_ptr_tran->create(atom_mu[J],atom_mu[I]);
                // cout << "I J: " << I << "  " << J << "   mu,nu: " << atom_mu[I] << "  " << atom_mu[J] << endl;
                for (size_t i_mu = 0; i_mu != basis_aux[I]; i_mu++)
                {
                    for (size_t i_nu = 0; i_nu != basis_aux[J]; i_nu++)
                    {
                        //(*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(J, i_nu), atom_mu_loc2glo(I, i_mu)); ////for aims
                        re(i_mu, i_nu) = vf_p.second(basis_aux.get_global_index(I, i_mu), basis_aux.get_global_index(J, i_nu)).real(); // for abacus
                        im(i_mu, i_nu) = vf_p.second(basis_aux.get_global_index(I, i_mu), basis_aux.get_global_index(J, i_nu)).imag();
                    }
                }

                if (is_cut_coulomb)
                {
                    driver::h.set_aux_cut_coulomb_k_atom_pair(iq, I, J, basis_aux[I], basis_aux[J], re.c, im.c, driver::opts.vq_threshold);
                }
                else
                {
                    driver::h.set_aux_bare_coulomb_k_atom_pair(iq, I, J, basis_aux[I], basis_aux[J], re.c, im.c, driver::opts.vq_threshold);
                }
                // if (I == J)
                // {
                //     (*vq_ptr).set_as_identity_matrix();
                // }

                // if ((*vq_ptr).real().absmax() >= threshold)
                // {
                //     coulomb_mat[I][J][qvec] = vq_ptr;
                //     vq_save++;
                // }
                // else
                // {
                //     vq_discard++;
                // }
            }
        }
    }
    profiler.stop("set_aux_coulomb_k_atom_pair_out");
    closedir(dir);
    dir = NULL;
    // cout << "vq threshold: " << threshold << endl;
    // cout << "vq_save:    " << vq_save << endl;
    // cout << "vq_dicard:  " << vq_discard << endl;
    // cout << "  Vq_dim   " << coulomb_mat.size() << "    " << coulomb_mat[0].size() << "   " << coulomb_mat[0][0].size() << endl;
    // for (auto &irk : irk_weight)
    // {
    //     cout << " irk_vec and weight: " << irk.first << "  " << irk.second << endl;
    // }
    // cout << "Finish read aims vq" << endl;
    return vq_discard;
}


static int handle_Vq_row_file(const string &file_path, double threshold,
        librpa_int::atom_mapping<std::map<int, std::shared_ptr<librpa_int::ComplexMatrix>>>::pair_t_old &coulomb,
        const std::vector<atpair_t> &local_atpair, bool binary)
{
    using librpa_int::ComplexMatrix;
    // cout << "Begin to read aims vq_real from " << file_path << endl;
    ifstream infile;
    int n_irk_points_local;
    int n_irk_points;

    auto ds = librpa_int::api::get_dataset_instance(driver::h.get_c_handler());
    const auto &basis_aux = ds->basis_aux;
    const auto atom_mu_part_range = basis_aux.get_part_range();

    if (binary)
    {
        infile.open(file_path, std::ios::in | std::ios::binary);
        infile.read((char *) &n_irk_points, sizeof(int));
        infile.read((char *) &n_irk_points_local, sizeof(int));
    }
    else
    {
        infile.open(file_path);
        infile >> n_irk_points;
    }
    if (!infile.good()) return 1;

    if (binary)
    {
        std::set<int> coulomb_row_need;
        for (const auto &[I, _]: local_atpair)
        {
            const auto brow = atom_mu_part_range[I];
            const auto nb = basis_aux[I];
            for (size_t ir = 0; ir < nb; ir++)
            {
                coulomb_row_need.insert(brow + ir);
            }
        }

        int nbasbas, brow, erow, bcol, ecol, iq;
        double q_weight;
        for (int i_irk = 0; i_irk < n_irk_points_local; i_irk++)
        {
            infile.read((char *) &nbasbas, sizeof(int));
            infile.read((char *) &brow, sizeof(int));
            infile.read((char *) &erow, sizeof(int));
            infile.read((char *) &bcol, sizeof(int));
            infile.read((char *) &ecol, sizeof(int));
            infile.read((char *) &iq, sizeof(int));
            infile.read((char *) &q_weight, sizeof(double));

            brow--;
            erow--;
            bcol--;
            ecol--;
            iq--;

            for (const auto &ap : local_atpair)
            {
                auto I = ap.first;
                auto J = ap.second;
                if (coulomb[I][J].count(iq) == 0)
                {
                    std::shared_ptr<ComplexMatrix> vq_ptr = std::make_shared<ComplexMatrix>();
                    vq_ptr->create(basis_aux[I], basis_aux[J]);
                    coulomb[I][J][iq] = vq_ptr;
                }
            }

            const auto ncol = ecol - bcol + 1;

            for (int i_mu = brow; i_mu <= erow; i_mu++)
            {
                std::vector<std::complex<double>> tmp_row(ncol);
                infile.read((char *) tmp_row.data(), 2 * ncol * sizeof(double));

                if (coulomb_row_need.count(i_mu))
                {
                    int I_loc, mu_loc;
                    basis_aux.get_local_index(i_mu, I_loc, mu_loc);
                    for (auto &Jp : coulomb[I_loc])
                    {
                        int J = Jp.first;
                        int Jb = atom_mu_part_range[J];
                        int Je = atom_mu_part_range[J] + basis_aux[J] - 1;

                        if (ecol >= Jb && bcol < Je)
                        {
                            int start_point = (bcol <= Jb ? Jb : bcol);
                            int end_point = (ecol <= Je ? ecol : Je);
                            for (int i = start_point; i <= end_point; i++)
                            {
                                int J_loc, nu_loc;
                                basis_aux.get_local_index(i, J_loc, nu_loc);
                                // printf("|i: %d   J: %d   J_loc: %d, nu_loc:
                                // %d\n",i,J,J_loc,nu_loc);
                                assert(J == J_loc);
                                (*coulomb[I_loc][J_loc][iq])(mu_loc, nu_loc) = tmp_row[i - bcol];
                            }
                        }
                    }
                }
            }
        }

    }
    else
    {
        string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3, vq_r, vq_i, q_num, q_weight;
        while (infile.peek() != EOF)
        {
            infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col;
            if (infile.peek() == EOF)
                break;
            if (!infile.good()) return 2;
            // cout << "vq range: " << begin_row << " ~ " << end_row << "  ,   " << begin_col << " ~ " << end_col << endl;
            infile >> q_num >> q_weight;
            if (!infile.good()) return 3;
            // int mu = stoi(nbasbas);
            // int nu = stoi(nbasbas);
            int brow = stoi(begin_row) - 1;
            int erow = stoi(end_row) - 1;
            int bcol = stoi(begin_col) - 1;
            int ecol = stoi(end_col) - 1;
            int iq = stoi(q_num) - 1;
            //cout<<file_path<<" iq:"<<iq<<"  qweight:"<<stod(q_weight)<<endl;

            //skip empty coulumb_file
            if((erow-brow<=0) || (ecol-bcol<=0) || iq<0)
                return 4;

            for(const auto &ap:local_atpair)
            {
                auto I=ap.first;
                auto J=ap.second;
                if(!coulomb[I][J].count(iq))
                {
                    std::shared_ptr<ComplexMatrix> vq_ptr = std::make_shared<ComplexMatrix>();
                    vq_ptr->create(basis_aux[I], basis_aux[J]);
                    // cout<<"  create  IJ: "<<I<<"  "<<J<<"   "<<atom_mu[I]<<"  "<<atom_mu[J];
                    coulomb[I][J][iq]=vq_ptr;
                }
            }   

            std::set<int> coulomb_row_need;
            for (auto &[I, _] : coulomb)
            {
                const int st = atom_mu_part_range[I];
                const int ed = atom_mu_part_range[I] + basis_aux[I];
                for (int ir = st; ir < ed; ir++) coulomb_row_need.insert(ir);
            }

            //printf("   |process %d, coulomb_begin:  %d, size: %d\n",para_mpi.get_myid(),*coulomb_row_need.begin(),coulomb_row_need.size());
            for (int i_mu = brow; i_mu <= erow; i_mu++)
            {
                std::vector<std::complex<double>> tmp_row(ecol-bcol+1);
                for (int i_nu = bcol; i_nu <= ecol; i_nu++)
                {
                    infile >> vq_r >> vq_i;
                    if (!infile.good()) return 4;

                    tmp_row[i_nu-bcol] = std::complex<double>(stod(vq_r), stod(vq_i)); // for abacus

                }
                if(coulomb_row_need.count(i_mu))
                {
                    int I_loc,mu_loc;
                    basis_aux.get_local_index(i_mu, I_loc,mu_loc);
                    // int bI=atom_mu_part_range[I_loc];
                    for(auto &Jp:coulomb[I_loc] )
                    {
                        int J=Jp.first;
                        int Jb=atom_mu_part_range[J];
                        int Je=atom_mu_part_range[J]+basis_aux[J]-1;

                        if(ecol>=Jb && bcol<Je)
                        {
                            int start_point = ( bcol<=Jb ? Jb:bcol);
                            int end_point = (ecol<=Je? ecol:Je);
                            for(int i=start_point;i<=end_point;i++)
                            {
                                int J_loc, nu_loc;
                                basis_aux.get_local_index(i,J_loc, nu_loc);
                                //printf("|i: %d   J: %d   J_loc: %d, nu_loc: %d\n",i,J,J_loc,nu_loc);
                                assert(J==J_loc);
                                (*coulomb[I_loc][J_loc][iq])(mu_loc,nu_loc)=tmp_row[i-bcol];
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}


size_t read_Vq_row(const string &dir_path, const string &vq_fprefix, double threshold,
        const std::vector<atpair_t> &local_atpair, bool is_cut_coulomb)
{
    using std::cout;
    using std::endl;
    using librpa_int::ComplexMatrix;
    using librpa_int::atom_mapping;
    using namespace librpa_int::global;

    cout << "Begin READ_Vq_Row" << endl;
    std::set<int> local_I_set;
    for(auto &lap:local_atpair)
    {
        local_I_set.insert(lap.first);
        local_I_set.insert(lap.second);
    }

    size_t vq_save = 0;
    size_t vq_discard = 0;
    atom_mapping<std::map<int, std::shared_ptr<ComplexMatrix>>>::pair_t_old coulomb;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    std::vector<string> files;
    bool binary;
    bool binary_checked = false;

    //map<Vector3_Order<double>, ComplexMatrix> Vq_full;
    profiler.start("handle_Vq_row_file");
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find(vq_fprefix) == 0)
        {
            string file_path = dir_path + fm;
            if (!binary_checked)
            {
                binary = check_coulomb_file_binary(file_path);
                binary_checked = true;
                if (librpa_int::global::myid_global == 0)
                {
                    if (binary)
                    {
                        cout << "Unformatted binary V files detected" << endl;
                    }
                    else
                    {
                        cout << "ASCII format V files detected" << endl;
                    }
                }
            }
            handle_Vq_row_file(file_path, threshold, coulomb, local_atpair, binary);
        }
    }
    profiler.stop("handle_Vq_row_file");

    // MYZ: now the map coulomb contains the complete atom-pair matrix.
    // Call the API to parse the data.
    // To reduce memory consumption during this process, we erase the data in temporary object once it is parsed.
    auto it_I = coulomb.begin();
    profiler.start("set_aux_coulomb_k_atom_pair_out");
    while (it_I != coulomb.end())
    {
        auto I = it_I->first;
        auto it_J = it_I->second.begin();
        while (it_J != it_I->second.end())
        {
            auto J = it_J->first;
            auto it_iq = it_J->second.begin();
            while (it_iq != it_J->second.end())
            {
                auto iq = it_iq->first;
                auto &vq_ptr = it_iq->second;
                if (is_cut_coulomb)
                {
                    driver::h.set_aux_cut_coulomb_k_atom_pair(iq, I, J, vq_ptr->nr, vq_ptr->nc, vq_ptr->real().c, vq_ptr->imag().c, threshold);
                }
                else
                {
                    driver::h.set_aux_bare_coulomb_k_atom_pair(iq, I, J, vq_ptr->nr, vq_ptr->nc, vq_ptr->real().c, vq_ptr->imag().c, threshold);
                }
                it_iq = it_J->second.erase(it_iq);
            }
            it_J = it_I->second.erase(it_J);
        }
        it_I = coulomb.erase(it_I);
    }
    profiler.stop("set_aux_coulomb_k_atom_pair_out");

    // cout << "FINISH coulomb files reading!" << endl;

    closedir(dir);
    dir = NULL;

    // ofstream fs;
    // std::stringstream ss;
    // ss<<"out_coulomb_rank_"<<para_mpi.get_myid()<<".txt";
    // fs.open(ss.str());
    // for(auto &Ip:coulomb_mat)
    // {
    //     for(auto &Jp:Ip.second)
    //         for(auto &qp:Jp.second)
    //         {
    //             std::stringstream sm;
    //             sm<<"I,J "<<Ip.first<<"  "<<Jp.first;
    //             //printf("|process %d  I J: %d, %d\n",para_mpi.get_myid(), Ip.first,Jp.first);
    //             print_complex_matrix_file(sm.str().c_str(),(*qp.second),fs,false);
    //         }

    // }
    // fs.close();
    return vq_discard;
}


void erase_Cs_from_local_atp(atpair_R_mat_t &Cs, std::vector<atpair_t> &local_atpair)
{
    using namespace std;
    using namespace librpa_int;
    //erase no need Cs

    set<size_t> loc_atp_index;
    for(auto &lap:local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    std::vector<atom_t> Cs_first;
    for (const auto &Ip: Cs)
        Cs_first.push_back(Ip.first);
    for (const auto &I: Cs_first)
    {
        if(!loc_atp_index.count(I))
            Cs.erase(I);
    }
    // for(auto &Ip:Cs)
    //     if(!loc_atp_index.count(Ip.first))
    //     {
    //         Cs.erase(Ip.first);
    //     }
    release_free_mem();
    global::lib_printf("| process %d, size of Cs after erase: %lu\n", librpa_int::global::mpi_comm_global_h.myid, Cs.size());
}

void read_stru(const std::string &file_path)
{
    using namespace librpa_int;
    global::lib_printf_root("Reading structure file: %s\n", file_path.c_str());

    ifstream infile;
    infile.open(file_path);
    string x, y, z, tmp;

    std::vector<double> lat_mat(9);
    std::vector<double> G_mat(9);

    for (int i = 0; i < 3; i++)
    {
        infile >> x >> y >> z;
        lat_mat[i * 3] = stod(x);
        lat_mat[i * 3 + 1] = stod(y);
        lat_mat[i * 3 + 2] = stod(z);
    }

    for (int i = 0; i < 3; i++)
    {
        infile >> x >> y >> z;
        G_mat[i * 3] = stod(x);
        G_mat[i * 3 + 1] = stod(y);
        G_mat[i * 3 + 2] = stod(z);
    }

    driver::h.set_latvec_and_G(lat_mat.data(), G_mat.data());

    // Read coordinates of atoms
    infile >> driver::n_atoms;
    const auto n_atoms = driver::n_atoms;
    driver::atom_types.resize(n_atoms);
    std::vector<double> coords(n_atoms * 3);
    int type;
    for (size_t iat = 0; iat < n_atoms; iat++)
    {
        for (int i = 0; i < 3; i++) infile >> coords[3 * iat + i];
        infile >> type;
        driver::atom_types[iat] = type - 1;
    }
    // Parsed after lattice is set, so that the fractional coordinates are calculated
    driver::h.set_atoms(driver::atom_types, coords);

    // // Internal check
    // const auto ds = api::get_dataset_instance(driver::h.get_c_handler());
    // const auto &pbc = ds->pbc;
    // Matrix3 latG = pbc.latvec * pbc.G.Transpose();
    // cout << " lat * G^T" << endl;
    // latG.print(5);
}

void read_bz_sampling(const std::string &file_path)
{
    using namespace librpa_int;

    global::lib_printf_root("Reading Brillouin zone sampling file: %s\n", file_path.c_str());

    ifstream infile;
    infile.open(file_path);

    string x, y, z, tmp;

    int nk[3];
    for (int i = 0; i < 3; i++)
    {
        infile >> nk[i];
    }
    int nk_full, nk_ibz;
    infile >> nk_full >> nk_ibz;
    assert(nk_full == nk[0] * nk[1] * nk[2]);

    std::vector<double> kvecs(3 * nk_full);
    std::vector<int> map_ibzk(nk_full, -1);

    // kvec_c = new Vector3<double>[n_kpoints];
    for (int i = 0; i != nk_full; i++)
    {
        // id weight kfrac[3] kcart[3] ik_ibz map_ibz
        infile >> x >> y;
        infile >> x >> y >> z;
        infile >> kvecs[3 * i] >> kvecs[3 * i + 1] >> kvecs[3 * i + 2];
        infile >> tmp >> map_ibzk[i];
        map_ibzk[i] -= 1;
        // Save a copy in the driver for information printing
        Vector3_Order<double> kvec(kvecs[3 * i], kvecs[3 * i + 1], kvecs[3 * i + 2]);
        auto it = std::find(driver::ibz_kpoints.cbegin(), driver::ibz_kpoints.cend(), kvec);
        if (it == driver::ibz_kpoints.cend()) driver::ibz_kpoints.emplace_back(kvec);
    }
    infile.close();

    driver::n_ibz_kpoints = driver::ibz_kpoints.size();

    driver::h.set_kgrids_kvec(nk[0], nk[1], nk[2], kvecs.data());
    driver::h.set_ibz_mapping(map_ibzk);
}

void read_basis(const std::string &file_path)
{
    using namespace librpa_int;

    global::lib_printf_root("Reading basis information file: %s\n", file_path.c_str());
    ifstream infile;
    infile.open(file_path);

    int n_atoms = driver::atom_types.size();
    if (as_size(n_atoms) != driver::n_atoms)
        throw std::runtime_error("Number of atoms not consistent with the geometry file!");
    std::map<int, size_t> map_at_wfc;
    std::map<int, size_t> map_at_aux;
    std::vector<size_t> nbs_wfc(n_atoms);
    std::vector<size_t> nbs_aux(n_atoms);

    int ntypes, type;
    size_t n_wfc, n_aux;
    string kind_str;

    infile >> ntypes;
    // total basis, not used here
    infile >> n_wfc >> n_aux >> kind_str;

    for (int itype = 0; itype < ntypes; itype++)
    {
        infile >> type >> n_wfc >> n_aux;
        type--;
        map_at_wfc[type] = n_wfc;
        map_at_aux[type] = n_aux;
    }
    for (int iat = 0; iat < n_atoms; iat++)
    {
        auto type = driver::atom_types[iat];
        nbs_wfc[iat] = map_at_wfc.at(type);
        nbs_aux[iat] = map_at_aux.at(type);
    }

    // std::cout << "nbs_wfc " << nbs_wfc << std::endl;
    // std::cout << "nbs_aux " << nbs_aux << std::endl;

    driver::h.set_ao_basis_wfc(nbs_wfc);
    driver::h.set_ao_basis_aux(nbs_aux);

    infile.close();
}


std::vector<Vector3_Order<double>> read_band_kpath_info(const string &file_path, int &n_basis, int &n_states, int &n_spin)
{
    std::vector<Vector3_Order<double>> kfrac_band;

    ifstream infile;
    infile.open(file_path);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    string x, y, z;
    int n_kpoints_band;

    // Read dimensions in the first row
    infile >> x;
    n_basis = stoi(x);
    infile >> x;
    n_states = stoi(x);
    infile >> x;
    n_spin = stoi(x);
    infile >> x;
    n_kpoints_band = stoi(x);

    for (int i = 0; i < n_kpoints_band; i++)
    {
        infile >> x >> y >> z;
        kfrac_band.push_back({stod(x), stod(y), stod(z)});
    }

    infile.close();

    return kfrac_band;
}

MeanField read_meanfield_band(const string &dir_path, int n_basis, int n_states, int n_spin, int n_kpoints_band)
{
    MeanField mf_band(n_spin, n_kpoints_band, n_states, n_basis);
    std::string s1, s2, s3, s4, s5;

    for (int ik = 0; ik < n_kpoints_band; ik++)
    {
        // Load occupation weights and eigenvalues
        std::stringstream ss;
        ss << dir_path << "band_KS_eigenvalue_k_" << std::setfill('0') << std::setw(5) << ik + 1 << ".txt";
        ifstream infile;
        infile.open(ss.str());

        for (int i_spin = 0; i_spin < n_spin; i_spin++)
        {
            for (int i_state = 0; i_state < n_states; i_state++)
            {
                infile >> s1 >> s2 >> s3 >> s4 >> s5;
                mf_band.get_weight()[i_spin](ik, i_state) = stod(s3);
                mf_band.get_eigenvals()[i_spin](ik, i_state) = stod(s4);
            }
        }

        infile.close();

        // Load eigenvectors
        ss.str("");
        ss.clear();
        ss << dir_path << "band_KS_eigenvector_k_" << std::setfill('0') << std::setw(5) << ik + 1 << ".txt";
        infile.open(ss.str(), std::ios::in | std::ios::binary);

        for (int i_spin = 0; i_spin < n_spin; i_spin++)
        {
            auto &wfc = mf_band.get_eigenvectors()[i_spin][ik];
            wfc.create(n_states, n_basis);
            const size_t nbytes = n_basis * n_states * sizeof(std::complex<double>);
            infile.read((char *) mf_band.get_eigenvectors()[i_spin][ik].c, nbytes);
        }

        infile.close();
    }

    // TODO: Fermi energy is not set

    return mf_band;
}

std::vector<matrix> read_vxc_band(const string &dir_path, int n_states, int n_spin, int n_kpoints_band)
{
    std::vector<matrix> vxc_band(n_spin);
    for (int i_spin = 0; i_spin < n_spin; i_spin++)
    {
        vxc_band[i_spin].create(n_kpoints_band, n_states);
    }
    std::string s1, s2, s3;

    for (int ik = 0; ik < n_kpoints_band; ik++)
    {
        // Load occupation weights and eigenvalues
        std::stringstream ss;
        ss << dir_path << "band_vxc_k_" << std::setfill('0') << std::setw(5) << ik + 1 << ".txt";
        ifstream infile;
        infile.open(ss.str());
        ss.clear();

        for (int i_spin = 0; i_spin < n_spin; i_spin++)
        {
            for (int i_state = 0; i_state < n_states; i_state++)
            {
                infile >> s1 >> s2 >> s3;
                vxc_band[i_spin](ik, i_state) = stod(s3);
            }
        }

        infile.close();
    }
    return vxc_band;
}

void read_elsi_csc(const string &file_path, bool save_row_major, std::vector<double> &mat, int &n_basis, bool &is_real)
{
    ifstream infile;
    infile.open(file_path, std::ios::binary);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    // Read the whole buffer
    infile.seekg(0, std::ios::end);
    std::streampos size = infile.tellg();
    infile.seekg(0, std::ios::beg);
    std::vector<char> buffer(size);
    infile.read(buffer.data(), size);
    infile.close();

    int64_t header[16];
    std::memcpy(header, buffer.data(), 128);

    n_basis = header[3];
    int64_t nnz = header[5];
    // cout << n_basis << " " << nnz << endl;

    int64_t* col_ptr_raw = reinterpret_cast<int64_t*>(buffer.data() + 128);
    std::vector<int> col_ptr;
    col_ptr.assign(col_ptr_raw, col_ptr_raw + n_basis);
    // Trailing column index to mark the end. +1 for index starting from 1 in ELSI CSC
    col_ptr.push_back(nnz + 1);

    int32_t* row_idx_raw = reinterpret_cast<int32_t*>(buffer.data() + 128 + n_basis * 8);

    char* nnz_val_raw = buffer.data() + 128 + n_basis * 8 + nnz * 4;
    double* nnz_val_double = reinterpret_cast<double*>(nnz_val_raw);

    if (header[2] == 0)
    {
        // Real valued
        is_real = true;
        mat.resize(n_basis * n_basis);
    } else {
        // Complex valued
        is_real = false;
        mat.resize(2 * n_basis * n_basis);
    }

    for (auto col = 0; col < n_basis; ++col) {
        for (auto idx = col_ptr[col]; idx < col_ptr[col + 1]; ++idx) {
            int row = row_idx_raw[idx - 1] - 1;
            int index = save_row_major ? row * n_basis + col : col * n_basis + row;
            // cout << idx - 1 << " " << col << " " << row << " " << index << endl;
            if (is_real)
            {
                mat[index] = nnz_val_double[idx - 1];
            }
            else
            {
                mat[2 * index] = nnz_val_double[2 * idx - 2];
                mat[2 * index + 1] = nnz_val_double[2 * idx - 1];
            }
        }
    }
}
