#include "read_aims.h"
#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include "atoms.h"
#include "ri.h"
#include "input.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
/* using std::stod; */

void READ_AIMS_BAND(const string &file_path, MeanField &mf)
{
    cout << "Begin to read aims-band_out" << endl;
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

    cout << "| number of spins: " << mf.get_n_spins() << endl
         << "| number of k-points: " << mf.get_n_kpoints() << endl
         << "| number of bands: " << mf.get_n_bands() << endl
         << "| number of NAOs: " << mf.get_n_aos() << endl;

    auto & eskb = mf.get_eigenvals();
    auto & wg = mf.get_weight();

    for (int is = 0; is != n_spins; is++)
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            infile >> ks >> ss;
            // cout<<ik<<is<<endl;
            int k_index = stoi(ks) - 1;
            int s_index = stoi(ss) - 1;
            for (int i = 0; i != n_bands; i++)
            {
                infile >> a >> ws >> es >> d;
                // cout<<a<<b<<c<<d<<endl;
                wg[is](k_index, i) = stod(ws) / n_kpoints; // different with abacus!
                eskb[is](k_index, i) = stod(es) * 2;
            }
        }
    cout << "efermi (Ha):" << efermi << endl;
    cout << "Success read aims-band_out" << endl;
    double emin, emax, gap;
    gap = mf.get_E_min_max(emin, emax);
    cout << "| Minimal transition energy (Ha): " << emin << endl
         << "| Maximal transition energy (Ha): " << emax << endl
         << "| Band gap (Ha): " << gap << endl;
}

void READ_AIMS_EIGENVECTOR(const string &dir_path, MeanField &mf)
{
    // cout<<"Begin to read aims eigenvecor"<<endl;
    assert(mf.get_n_spins() == 1);

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
    cout << "Finish read KS_eignvector! " << endl;
}


void handle_KS_file(const string &file_path, MeanField &mf)
{
    // cout<<file_path<<endl;
    ifstream infile;
    cout << "Reading eigenvector from file " << file_path << endl;
    infile.open(file_path);
    int ik;
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
    }}

void READ_AIMS_Cs(const string &dir_path, double threshold)
{
    cout << "Begin to read Cs" << endl;
    cout << "cs_threshold:  " << threshold << endl;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("Cs_data") == 0)
            handle_Cs_file(fm, threshold);
    }
    closedir(dir);
    dir = NULL;
    
    atom_mu_part_range.resize(atom_mu.size());
    atom_mu_part_range[0]=0;
    for(int I=1;I!=atom_mu.size();I++)
        atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    for(int i=0;i!=atom_mu_part_range.size();i++)
        cout<<" atom_mu_part_range ,i: "<<i<<"    "<<atom_mu_part_range[i]<<endl;

    cout << "Finish read Cs" << endl;
}

void handle_Cs_file(const string &file_path, double threshold)
{
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
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

        // if (box == Vector3_Order<int>({0, 0, 1}))continue;
        if ((*cs_ptr).absmax() >= threshold)
            Cs[ia1][ia2][box] = cs_ptr;

        // cout<<" READ Cs, INDEX:  "<<ia1<<"   "<<ia2<<"   "<<box<<"   "<<(*Cs.at(ia1).at(ia2).at(box))(n_i*n_j-1,n_mu-1)<<endl;
    }
}

void READ_AIMS_Vq(const string &dir_path, double threshold)
{
    int vq_save = 0;
    int vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    map<Vector3_Order<double>, ComplexMatrix> Vq_full;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("coulomb_mat") == 0)
        {
            handle_Vq_file(fm, threshold, Vq_full);
        }
    }
    cout << "FINISH coulomb files reading!" << endl;
    for (auto &vf_p : Vq_full)
    {
        auto qvec = vf_p.first;
        cout << "Qvec:" << qvec << endl;
        for (int I = 0; I != atom_mu.size(); I++)
            for (int J = 0; J != atom_mu.size(); J++)
            {
                if (I > J)
                    continue;
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);
                // vq_ptr_tran->create(atom_mu[J],atom_mu[I]);
                cout << "I J: " << I << "  " << J << "   mu,nu: " << atom_mu[I] << "  " << atom_mu[J] << endl;
                for (int i_mu = 0; i_mu != atom_mu[I]; i_mu++)
                {

                    for (int i_nu = 0; i_nu != atom_mu[J]; i_nu++)
                    {
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(J, i_nu), atom_mu_loc2glo(I, i_mu)); ////for aims
                        //(*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                // if (I == J)
                // {
                //     (*vq_ptr).set_as_identity_matrix();
                // }

                if ((*vq_ptr).real().absmax() >= threshold)
                {
                    Vq[I][J][qvec] = vq_ptr;
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
    cout << "vq threshold: " << threshold << endl;
    cout << "vq_save:    " << vq_save << endl;
    cout << "vq_dicard:  " << vq_discard << endl;
    cout << "  Vq_dim   " << Vq.size() << "    " << Vq[0].size() << "   " << Vq[0][0].size() << endl;
    for (auto &irk : irk_weight)
    {
        cout << " irk_vec and weight: " << irk.first << "  " << irk.second << endl;
        // Cal_Periodic_Chi0::print_complex_matrix("full_Vq",Vq_full.at(irk.first));
    }
    cout << "Finish read aims vq" << endl;
}

void handle_Vq_file(const string &file_path, double threshold, map<Vector3_Order<double>, ComplexMatrix> &Vq_full)
{
    cout << "Begin to read aims vq_real from " << file_path << endl;
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

    while (infile.peek() != EOF)
    {
        infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col;
        if (infile.peek() == EOF)
            break;
        // cout << "vq range: " << begin_row << " ~ " << end_row << "  ,   " << begin_col << " ~ " << end_col << endl;
        infile >> q_num >> q_weight;
        int mu = stoi(nbasbas);
        int nu = stoi(nbasbas);
        int brow = stoi(begin_row) - 1;
        int erow = stoi(end_row) - 1;
        int bcol = stoi(begin_col) - 1;
        int ecol = stoi(end_col) - 1;
        int iq = stoi(q_num) - 1;
        Vector3_Order<double> qvec(kvec_c[iq]);
        irk_weight.insert(pair<Vector3_Order<double>, double>(qvec, stod(q_weight)));
        if (!Vq_full.count(qvec))
        {
            Vq_full[qvec].create(mu, nu);
        }
        for (int i_mu = brow; i_mu <= erow; i_mu++)
            for (int i_nu = bcol; i_nu <= ecol; i_nu++)
            {
                infile >> vq_r >> vq_i;
                Vq_full[qvec](i_nu, i_mu) = complex<double>(stod(vq_r), stod(vq_i)); // for FHI-aims
                // Vq_full[qvec](i_mu, i_nu) = complex<double>(stod(vq_r), stod(vq_i)); // for abacus
            }
    }
}

// TODO: implement the wrapper of all input readers
void read_aims(MeanField &mf)
{
}
