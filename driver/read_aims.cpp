#include "read_aims.h"
// #include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <algorithm>
#include <unordered_map>
#if defined(__MACH__)
#include <malloc/malloc.h> // for malloc_zone_pressure_relief and malloc_default_zone
#else
#include <malloc.h>
#endif
#include "atoms.h"
#include "atomic_basis.h"
#include "ri.h"
#include "pbc.h"
#include "constants.h"
#include "envs_mpi.h"
#include "envs_io.h"
#include "stl_io_helper.h"

// using std::cout;
// using std::endl;
using std::ifstream;
using std::string;
/* using std::stod; */

void READ_AIMS_BAND(const string &file_path, MeanField &mf)
{
    // cout << "Begin to read aims-band_out" << endl;
    ifstream infile;
    infile.open(file_path);
    string ks, ss, a, ws, es, d;
    int n_kpoints, n_spins, n_bands, n_aos;
    double efermi;
    infile >> n_kpoints;
    infile >> n_spins;
    infile >> n_bands;
    infile >> n_aos;
    infile >> efermi;
    efermi *= 2;
    mf.get_efermi() = efermi;   
    // efermi=0.2;

    mf.set(n_spins, n_kpoints, n_bands, n_aos);

    auto & eskb = mf.get_eigenvals();
    auto & wg = mf.get_weight();
    //cout<<"|eskb: "<<endl;
    for (int ik = 0; ik != n_kpoints; ik++)
    {
        for (int is = 0; is != n_spins; is++)
        {
            infile >> ks >> ss;
            //cout<<ik<<is<<endl;
            int k_index = stoi(ks) - 1;
            // int s_index = stoi(ss) - 1;
            for (int i = 0; i != n_bands; i++)
            {
                infile >> a >> ws >> es >> d;
                // cout<<a<<b<<c<<d<<endl;
                wg[is](k_index, i) = stod(ws) / n_kpoints; // different with abacus!
                eskb[is](k_index, i) = stod(es) * 2;
                //cout<<" i_band: "<<i<<"    eskb: "<<eskb[is](k_index, i)<<endl;
            }
        }
    }
    // for (int is = 0; is != n_spins; is++)
    //     print_matrix("eskb_mat",eskb[is]);
}

int read_vxc(const string &file_path, std::vector<matrix> &vxc)
{
    // cout << "Begin to read aims-band_out" << endl;
    ifstream infile;
    infile.open(file_path);
    double ha, ev;
    int n_spins, n_kpoints, n_states;
    int retcode;

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

void READ_AIMS_EIGENVECTOR(const string &dir_path, MeanField &mf)
{
    // cout<<"Begin to read aims eigenvecor"<<endl;
    //assert(mf.get_n_spins() == 1);

    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("KS_eigenvector") == 0)
        {
            handle_KS_file(fm, mf);
        }
    }
    closedir(dir);
    dir = NULL;
    //auto tmp_wfc=mf.get_eigenvectors();
    // for(int is=0;is!=mf.get_n_spins();is++)
    //     print_complex_matrix("wfc ",tmp_wfc.at(is).at(0));
    // cout << "Finish read KS_eignvector! " << endl;
}


void handle_KS_file(const string &file_path, MeanField &mf)
{
    // cout<<file_path<<endl;
    ifstream infile;
    // cout << "Reading eigenvector from file " << file_path << endl;
    infile.open(file_path);
    // int ik;
    string rvalue, ivalue, kstr;
    auto & wfc = mf.get_eigenvectors();
    while (infile.peek() != EOF)
    {
        infile >> kstr;
        // cout<<"     ik: "<<ik<<endl;
        if (infile.peek() == EOF)
            break;
        // for aims !!!
        for (int iw = 0; iw != mf.get_n_aos(); iw++)
            for (int ib = 0; ib != mf.get_n_bands(); ib++)
                for (int is = 0; is != mf.get_n_spins(); is++)
                {
                    // cout<<iw<<ib<<is<<ik;
                    infile >> rvalue >> ivalue;
                    // cout<<rvalue<<ivalue<<endl;
                    wfc.at(is).at(stoi(kstr) - 1)(ib, iw) = complex<double>(stod(rvalue), stod(ivalue));
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
}

size_t READ_AIMS_Cs(const string &dir_path, double threshold,const vector<atpair_t> &local_atpair, bool binary)
{
    size_t cs_discard = 0;
    // cout << "Begin to read Cs" << endl;
    // cout << "cs_threshold:  " << threshold << endl;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("Cs_data") == 0)
        {
            if (binary)
            {
                cs_discard += handle_Cs_file_binary(fm, threshold, local_atpair);
            }
            else
            {
                cs_discard += handle_Cs_file(fm, threshold, local_atpair);
            }
        }
    }
    closedir(dir);
    dir = NULL;
    // initialize basis set object
    LIBRPA::atomic_basis_wfc.set(atom_nw);
    LIBRPA::atomic_basis_abf.set(atom_mu);
    
    // atom_mu_part_range.resize(atom_mu.size());
    // atom_mu_part_range[0]=0;
    // for(int I=1;I!=atom_mu.size();I++)
    //     atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    // N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    init_N_all_mu();

    // for(int i=0;i!=atom_mu_part_range.size();i++)
    //     cout<<" atom_mu_part_range ,i: "<<i<<"    "<<atom_mu_part_range[i]<<endl;

    // cout << "Finish read Cs" << endl;
    return cs_discard;
}

size_t READ_AIMS_Cs_evenly_distribute(const string &dir_path, double threshold, int myid, int nprocs, bool binary)
{
    size_t cs_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    unordered_map<string, vector<size_t>> files_Cs_ids_this_proc;
    int Cs_keep_total = 0;

    while ((ptr = readdir(dir)) != NULL)
    {
        string fn(ptr->d_name);
        if (fn.find("Cs_data") == 0)
        {
            files.push_back(fn);
            std::vector<size_t> ids_keep_this_file;
            if (binary)
            {
                ids_keep_this_file = handle_Cs_file_binary_dry(fn, threshold);
            }
            else
            {
                ids_keep_this_file = handle_Cs_file_dry(fn, threshold);
            }
            for (int id = 0; id < ids_keep_this_file.size(); id++)
            {
                int id_global = id + Cs_keep_total;
                if (id_global % nprocs == myid) files_Cs_ids_this_proc[fn].push_back(ids_keep_this_file[id]);
            }
            Cs_keep_total += ids_keep_this_file.size();
        }
    }
    closedir(dir);
    dir = NULL;
    if (myid == 0) printf("Finished Cs filtering\n");

    for (const auto& fn_ids: files_Cs_ids_this_proc)
    {
        LIBRPA::envs::ofs_myid << fn_ids.first << " " << fn_ids.second << endl;
        if (binary)
        {
            cs_discard += handle_Cs_file_binary_by_ids(fn_ids.first, threshold, fn_ids.second);
        }
        else
        {
            cs_discard += handle_Cs_file_by_ids(fn_ids.first, threshold, fn_ids.second);
        }
    }

    // initialize basis set object
    LIBRPA::atomic_basis_wfc.set(atom_nw);
    LIBRPA::atomic_basis_abf.set(atom_mu);
    
    atom_mu_part_range.resize(atom_mu.size());
    atom_mu_part_range[0]=0;
    for(int I=1;I!=atom_mu.size();I++)
        atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    return cs_discard;
}

void get_natom_ncell_from_first_Cs_file(int &n_atom, int &n_cell, const string &dir_path, bool binary)
{
    // cout<<file_path<<endl;
    ifstream infile;

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
            file_path = fn;
            break;
        }
    }
    if (file_path == "")
        throw std::runtime_error("Cs_data file is not found under dir_path: " + dir_path);

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
    while(ifs >> omega >> re >> im)
    {
        omegas.push_back(omega);
        dielec_func_imagfreq.push_back(re);
    }
    ifs.close();
}

size_t handle_Cs_file(const string &file_path, double threshold, const vector<atpair_t> &local_atpair)
{
    
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
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    natom = stoi(natom_s);
    ncell = stoi(ncell_s);
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
        int ic1 = stoi(ic_1);
        int ic2 = stoi(ic_2);
        int ic3 = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        atom_nw.insert(pair<int, int>(ia1, n_i));
        atom_mu.insert(pair<int, int>(ia1, n_mu));
        Vector3_Order<int> box(ic1, ic2, ic3);
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
        if (loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold )
            Cs[ia1][ia2][box] = cs_ptr;
        else
            cs_discard++;
    }
    return cs_discard;
}

size_t handle_Cs_file_binary(const string &file_path, double threshold, const vector<atpair_t> &local_atpair)
{
    
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

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        int ic1 = dims[2];
        int ic2 = dims[3];
        int ic3 = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        atom_nw.insert(pair<int, int>(ia1, n_i));
        atom_mu.insert(pair<int, int>(ia1, n_mu));
        Vector3_Order<int> box(ic1, ic2, ic3);
        // cout<< ia1<<ia2<<box<<endl;

        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(n_i * n_j, n_mu);
        // cout<<cs_ptr->nr<<cs_ptr->nc<<endl;

        infile.read((char *) cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
        if (loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold )
            Cs[ia1][ia2][box] = cs_ptr;
        else
            cs_discard++;
    }
    return cs_discard;
}

std::vector<size_t> handle_Cs_file_dry(const string &file_path, double threshold)
{
    std::vector<size_t> Cs_ids_keep;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    natom = stoi(natom_s);
    ncell = stoi(ncell_s);

    size_t id = 0;
    while (infile.peek() != EOF)
    {
        infile >> ia1_s;
        if (infile.peek() == EOF)
            break;
        infile >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s >> j_s >> mu_s;
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        double maxval = -1.0;
        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> Cs_ele;
                    maxval = std::max(maxval, abs(stod(Cs_ele)));
                }
        LIBRPA::envs::ofs_myid << id << " (" << ic_1 << "," << ic_2 << "," << ic_3 << ") " << maxval << " keep? " << (maxval >= threshold) << endl;
        if (maxval >= threshold)
            Cs_ids_keep.push_back(id);
        id++;
    }
    LIBRPA::envs::ofs_myid << file_path << ": " << Cs_ids_keep << endl;
    infile.close();
    return Cs_ids_keep;
}

std::vector<size_t> handle_Cs_file_binary_dry(const string &file_path, double threshold)
{
    std::vector<size_t> Cs_ids_keep;
    ifstream infile;
    int dims[8];
    int n_apcell_file;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));

    for (int i_file = 0; i_file < n_apcell_file; i_file++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        int ic1 = dims[2];
        int ic2 = dims[3];
        int ic3 = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        // cout<< ia1<<ia2<<box<<endl;
        double maxval = -1.0;
        double Cs_read;
        for (int i = 0; i != n_i; i++)
        {
            for (int j = 0; j != n_j; j++)
            {
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile.read((char *) &Cs_read, sizeof(double));
                    maxval = std::max(maxval, abs(Cs_read));
                }
            }
        }
        LIBRPA::envs::ofs_myid << i_file << " (" << ic1 << "," << ic2 << "," << ic3 << ") " << maxval << " keep? " << (maxval >= threshold) << endl;
        if (maxval >= threshold)
            Cs_ids_keep.push_back(i_file);
    }
    LIBRPA::envs::ofs_myid << file_path << ": " << Cs_ids_keep << endl;
    infile.close();
    return Cs_ids_keep;
}

size_t handle_Cs_file_by_ids(const string &file_path, double threshold, const vector<size_t> &ids)
{
    size_t cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    natom = stoi(natom_s);
    ncell = stoi(ncell_s);
    /* cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl; */
    // for(int loop=0;loop!=natom*natom*ncell;loop++)
    size_t id = 0;
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
        int ic1 = stoi(ic_1);
        int ic2 = stoi(ic_2);
        int ic3 = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        atom_nw.insert(pair<int, int>(ia1, n_i));
        atom_mu.insert(pair<int, int>(ia1, n_mu));
        Vector3_Order<int> box(ic1, ic2, ic3);

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
            Cs[ia1][ia2][box] = cs_ptr;
        }
        else
        {
            double maxval = -1.0;
            for (int i = 0; i != n_i; i++)
                for (int j = 0; j != n_j; j++)
                    for (int mu = 0; mu != n_mu; mu++)
                    {
                        infile >> Cs_ele;
                        maxval = std::max(maxval, abs(stod(Cs_ele)));
                    }
            if (maxval < threshold) cs_discard++;
        }
        id++;
    }
    infile.close();
    return cs_discard;
}

size_t handle_Cs_file_binary_by_ids(const string &file_path, double threshold, const vector<size_t> &ids)
{
    ifstream infile;
    int dims[8];
    int n_apcell_file;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));
    size_t cs_discard = 0;

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *) &dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        int ic1 = dims[2];
        int ic2 = dims[3];
        int ic3 = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        atom_nw.insert(pair<int, int>(ia1, n_i));
        atom_mu.insert(pair<int, int>(ia1, n_mu));
        Vector3_Order<int> box(ic1, ic2, ic3);

        if (std::find(ids.cbegin(), ids.cend(), static_cast<size_t>(i)) != ids.cend())
        {
            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);
            infile.read((char *) cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
            Cs[ia1][ia2][box] = cs_ptr;
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

size_t READ_Vq_Full(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_mat_t &coulomb_mat)
{
    size_t vq_save = 0;
    size_t vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    map<Vector3_Order<double>, ComplexMatrix> Vq_full;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find(vq_fprefix) == 0)
        {
            int retcode = handle_Vq_full_file(fm, threshold, Vq_full);
            if (retcode != 0)
            {
                printf("Error encountered when reading %s, return code %d", fm.c_str(), retcode);
            }
        }
    }
    // cout << "FINISH coulomb files reading!" << endl;
    for (auto &vf_p : Vq_full)
    {
        auto qvec = vf_p.first;
        // cout << "Qvec:" << qvec << endl;
        for (int I = 0; I != atom_mu.size(); I++)
            for (int J = 0; J != atom_mu.size(); J++)
            {
                if (I > J)
                    continue;
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);
                // vq_ptr_tran->create(atom_mu[J],atom_mu[I]);
                // cout << "I J: " << I << "  " << J << "   mu,nu: " << atom_mu[I] << "  " << atom_mu[J] << endl;
                for (int i_mu = 0; i_mu != atom_mu[I]; i_mu++)
                {

                    for (int i_nu = 0; i_nu != atom_mu[J]; i_nu++)
                    {
                        //(*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(J, i_nu), atom_mu_loc2glo(I, i_mu)); ////for aims
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                // if (I == J)
                // {
                //     (*vq_ptr).set_as_identity_matrix();
                // }

                if ((*vq_ptr).real().absmax() >= threshold)
                {
                    coulomb_mat[I][J][qvec] = vq_ptr;
                    vq_save++;
                }
                else
                {
                    vq_discard++;
                }
            }
    }
    closedir(dir);
    dir = NULL;
    // cout << "vq threshold: " << threshold << endl;
    // cout << "vq_save:    " << vq_save << endl;
    // cout << "vq_dicard:  " << vq_discard << endl;
    // cout << "  Vq_dim   " << coulomb_mat.size() << "    " << coulomb_mat[0].size() << "   " << coulomb_mat[0][0].size() << endl;
    // for (auto &irk : irk_weight)
    // {
    //     cout << " irk_vec and weight: " << irk.first << "  " << irk.second << endl;
    //     Cal_Periodic_Chi0::print_complex_matrix("full_Vq",Vq_full.at(irk.first));
    // }
    // cout << "Finish read aims vq" << endl;
    return vq_discard;
}

int handle_Vq_full_file(const string &file_path, double threshold, map<Vector3_Order<double>, ComplexMatrix> &Vq_full)
{
    // cout << "Begin to read aims vq_real from " << file_path << endl;
    ifstream infile;
    infile.open(file_path);
    string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3, vq_r, vq_i, q_num, q_weight;
    // int nline=0;
    // while(!infile.eof())
    // {
    //     nline++;
    // }
    // cout<<"  nline:  "<<nline<<endl;
    infile >> n_irk_points;
    if (!infile.good())
        return 1;

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
        if((erow-brow<=0) || (ecol-bcol<=0) || iq<0 || iq> klist.size())
            return 4;
        Vector3_Order<double> qvec(kvec_c[iq]);
        // skip duplicate insert of k weight, since 
        if (irk_weight.count(qvec) == 0)
            irk_weight.insert(pair<Vector3_Order<double>, double>(qvec, stod(q_weight)));
        if (!Vq_full.count(qvec))
        {
            Vq_full[qvec].create(mu, nu);
        }
        for (int i_mu = brow; i_mu <= erow; i_mu++)
            for (int i_nu = bcol; i_nu <= ecol; i_nu++)
            {
                infile >> vq_r >> vq_i;
                //Vq_full[qvec](i_nu, i_mu) = complex<double>(stod(vq_r), stod(vq_i)); // for FHI-aims
                Vq_full[qvec](i_mu, i_nu) = complex<double>(stod(vq_r), stod(vq_i)); // for abacus
            }
    }
    return 0;
}


size_t READ_Vq_Row(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_mat_t &coulomb_mat, const vector<atpair_t> &local_atpair)
{
    cout<<"Begin READ_Vq_Row"<<endl;
    set<int> local_I_set;
    for(auto &lap:local_atpair)
    {
        local_I_set.insert(lap.first);
        local_I_set.insert(lap.second);
    }

    size_t vq_save = 0;
    size_t vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    //map<Vector3_Order<double>, ComplexMatrix> Vq_full;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find(vq_fprefix) == 0)
        {
            //handle_Vq_full_file(fm, threshold, Vq_full);
            handle_Vq_row_file(fm,threshold,coulomb_mat,local_atpair);
        }
    }
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

int handle_Vq_row_file(const string &file_path, double threshold, atpair_k_cplx_mat_t &coulomb, const vector<atpair_t> &local_atpair)
{
    // cout << "Begin to read aims vq_real from " << file_path << endl;
    ifstream infile;
    infile.open(file_path);
    string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3, vq_r, vq_i, q_num, q_weight;
    infile >> n_irk_points;
    if (!infile.good()) return 1;

    while (infile.peek() != EOF)
    {
        infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col;
        if (infile.peek() == EOF)
            break;
        if (!infile.good()) return 2;
        // cout << "vq range: " << begin_row << " ~ " << end_row << "  ,   " << begin_col << " ~ " << end_col << endl;
        infile >> q_num >> q_weight;
        if (!infile.good()) return 3;
        int mu = stoi(nbasbas);
        int nu = stoi(nbasbas);
        int brow = stoi(begin_row) - 1;
        int erow = stoi(end_row) - 1;
        int bcol = stoi(begin_col) - 1;
        int ecol = stoi(end_col) - 1;
        int iq = stoi(q_num) - 1;
        //cout<<file_path<<" iq:"<<iq<<"  qweight:"<<stod(q_weight)<<endl;

        //skip empty coulumb_file
        if((erow-brow<=0) || (ecol-bcol<=0) || iq<0 || iq> klist.size())
            return 4;

        Vector3_Order<double> qvec(kvec_c[iq]);
        // skip duplicate insert of k weight, since 
        if (irk_weight.count(qvec) == 0)
            irk_weight.insert(pair<Vector3_Order<double>, double>(qvec, stod(q_weight)));

        for(const auto &ap:local_atpair)
        {
            auto I=ap.first;
            auto J=ap.second;
            if(!coulomb[I][J].count(qvec))
            {
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);
                // cout<<"  create  IJ: "<<I<<"  "<<J<<"   "<<atom_mu[I]<<"  "<<atom_mu[J];
                coulomb[I][J][qvec]=vq_ptr;
            }
        }   

        set<int> coulomb_row_need;
        for(auto &Ip:coulomb)
            for(int ir=atom_mu_part_range[Ip.first];ir!=atom_mu_part_range[Ip.first]+atom_mu[Ip.first];ir++)
                coulomb_row_need.insert(ir);

        //printf("   |process %d, coulomb_begin:  %d, size: %d\n",para_mpi.get_myid(),*coulomb_row_need.begin(),coulomb_row_need.size());
        for (int i_mu = brow; i_mu <= erow; i_mu++)
        {
            vector<complex<double>> tmp_row(ecol-bcol+1);
            for (int i_nu = bcol; i_nu <= ecol; i_nu++)
            {
                infile >> vq_r >> vq_i;
                if (!infile.good()) return 4;
                
                tmp_row[i_nu-bcol] = complex<double>(stod(vq_r), stod(vq_i)); // for abacus
                
            }
            if(coulomb_row_need.count(i_mu))
            {
                int I_loc,mu_loc;
                I_loc=atom_mu_glo2loc(i_mu,mu_loc);
                int bI=atom_mu_part_range[I_loc];
                for(auto &Jp:coulomb[I_loc] )
                {
                    auto J=Jp.first;
                    int Jb=atom_mu_part_range[J];
                    int Je=atom_mu_part_range[J]+atom_mu[J]-1;
                    
                    if(ecol>=Jb && bcol<Je)
                    {
                        int start_point = ( bcol<=Jb ? Jb:bcol);
                        int end_point = (ecol<=Je? ecol:Je);
                        for(int i=start_point;i<=end_point;i++)
                        {
                            int J_loc, nu_loc;
                            J_loc=atom_mu_glo2loc(i,nu_loc);
                            //printf("|i: %d   J: %d   J_loc: %d, nu_loc: %d\n",i,J,J_loc,nu_loc);
                            assert(J==J_loc);
                            (*coulomb[I_loc][J_loc][qvec])(mu_loc,nu_loc)=tmp_row[i-bcol];
                        }
                    }
                }
            }
        }
    }
    return 0;
}

void erase_Cs_from_local_atp(atpair_R_mat_t &Cs, vector<atpair_t> &local_atpair)
{
    //erase no need Cs
    
    set<size_t> loc_atp_index;
    for(auto &lap:local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    vector<atom_t> Cs_first;
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
    #ifndef __MACH__
    malloc_trim(0);
    #else
    malloc_zone_pressure_relief(malloc_default_zone(), 0);
    #endif
    LIBRPA::utils::lib_printf("| process %d, size of Cs after erase: %lu\n", LIBRPA::envs::mpi_comm_global_h.myid, Cs.size());
}

void READ_AIMS_STRU(const int& n_kpoints, const std::string &file_path)
{
    // cout << "Begin to read aims stru" << endl;
    ifstream infile;
    string x, y, z;
    infile.open(file_path);
    infile >> x >> y >> z;
    latvec.e11 = stod(x);
    latvec.e12 = stod(y);
    latvec.e13 = stod(z);
    infile >> x >> y >> z;
    latvec.e21 = stod(x);
    latvec.e22 = stod(y);
    latvec.e23 = stod(z);
    infile >> x >> y >> z;
    latvec.e31 = stod(x);
    latvec.e32 = stod(y);
    latvec.e33 = stod(z);
    latvec /= ANG2BOHR;
    // latvec.print();
    lat_array[0] = {latvec.e11,latvec.e12,latvec.e13};
    lat_array[1] = {latvec.e21,latvec.e22,latvec.e23};
    lat_array[2] = {latvec.e31,latvec.e32,latvec.e33};

    infile >> x >> y >> z;
    G.e11 = stod(x);
    G.e12 = stod(y);
    G.e13 = stod(z);
    infile >> x >> y >> z;
    G.e21 = stod(x);
    G.e22 = stod(y);
    G.e23 = stod(z);
    infile >> x >> y >> z;
    G.e31 = stod(x);
    G.e32 = stod(y);
    G.e33 = stod(z);

    G /= TWO_PI;
    G *= ANG2BOHR;
    // G.print();
    // Matrix3 latG = latvec * G.Transpose();
    // cout << " lat * G^T" << endl;
    // latG.print();
    infile >> x >> y >> z;
    kv_nmp[0] = stoi(x);
    kv_nmp[1] = stoi(y);
    kv_nmp[2] = stoi(z);
    kvec_c = new Vector3<double>[n_kpoints];
    for (int i = 0; i != n_kpoints; i++)
    {
        infile >> x >> y >> z;
        kvec_c[i] = {stod(x), stod(y), stod(z)};
        kvec_c[i] *= (ANG2BOHR / TWO_PI);
        // cout << "kvec [" << i << "]: " << kvec_c[i] << endl;
        Vector3_Order<double> k(kvec_c[i]);
        klist.push_back(k);
        kfrac_list.push_back(latvec * k);
    }
    for (int i = 0; i != n_kpoints; i++)
    {
        infile >> x;
        int id_irk = stoi(x) - 1;
        irk_point_id_mapping.push_back(id_irk);
        map_irk_ks[klist[id_irk]].push_back(klist[i]);
    }
}

// TODO: implement the wrapper of all input readers
void read_aims(MeanField &mf)
{
}
