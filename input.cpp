#include "input.h"
#include <fstream>
#include <iostream>
#include <dirent.h>
#include "cal_periodic_chi0.h"
#include "constants.h"
using namespace std;
double cs_threshold = 1E-6;
int NBANDS = 0;
int NLOCAL = 0;
int NSPIN = 0;
int natom = 0;
int ncell = 0;
matrix wg;
double **ekb;
double efermi = 0;
int n_kpoints = 0;
int n_irk_points = 0;
int kv_nmp[3] = {1, 1, 1};
Vector3<double> *kvec_c;
Matrix3 latvec;
Matrix3 G;
// vector<ComplexMatrix> wfc_k;
map<size_t, map<size_t, map<Vector3_Order<int>, std::shared_ptr<matrix>>>> Cs;
map<size_t, map<size_t, map<Vector3_Order<double>, std::shared_ptr<ComplexMatrix>>>> Vq;
map<Vector3_Order<double>, double> irk_weight;
map<int, int> atom_nw;
map<int, int> atom_mu;

void READ_AIMS_BAND(const std::string &file_path)
{
    cout << "Begin to read aims-band_out" << endl;
    ifstream infile;
    infile.open(file_path);
    string ik, is, a, b, c, d;
    infile >> n_kpoints;
    infile >> NSPIN;
    infile >> NBANDS;
    infile >> NLOCAL;
    infile >> efermi;
    efermi *= 2;
    // efermi=0.2;
    cout << n_kpoints << NSPIN << NLOCAL << endl;
    wg.create(n_kpoints, NLOCAL);
    ekb = new double *[n_kpoints];
    for (int ik = 0; ik != n_kpoints; ik++)
    {
        ekb[ik] = new double[NBANDS];
    }

    for (int ks = 0; ks != n_kpoints; ks++)
    {
        infile >> ik >> is;
        // cout<<ik<<is<<endl;
        int k_index = stoi(ik) - 1;
        int s_index = stoi(is) - 1;
        for (int i = 0; i != NBANDS; i++)
        {
            infile >> a >> b >> c >> d;
            // cout<<a<<b<<c<<d<<endl;
            wg(k_index, i) = stod(b) / n_kpoints; // different with abacus!
            ekb[k_index][i] = stod(c) * 2;
        }
    }
    cout << "efermi :" << efermi << endl;
    cout << "Success read aims-band_out" << endl;
}

void READ_AIMS_STRU(const std::string &file_path)
{
    cout << "Begin to read aims stru" << endl;
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
    latvec.print();

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
    G.print();
    Matrix3 latG = latvec * G;
    cout << " lat * G" << endl;
    latG.print();
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
        cout << kvec_c[i] << endl;
    }
}

void READ_AIMS_EIGENVECTOR(const std::string &file_path, vector<ComplexMatrix> &wfc_k)
{
    // cout<<"Begin to read aims eigenvecor"<<endl;
    assert(NSPIN == 1);
    // cout<<" n_kpoints:" <<n_kpoints<<endl;
    wfc_k.resize(n_kpoints);
    // cout<<"wfc_k size: "<<wfc_k.size()<<endl;
    for (int ik = 0; ik != n_kpoints; ik++)
        wfc_k[ik].create(NBANDS, NLOCAL);
    // cout<<"  create wfc_k"<<endl;

    struct dirent *ptr;
    DIR *dir;
    dir = opendir(file_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("KS_eigenvector") == 0)
        {
            handle_KS_file(fm, wfc_k);
        }
    }
    closedir(dir);
    dir = NULL;
    cout << "Finish read KS_eignvector! " << endl;
}

void handle_KS_file(const std::string &file_path, vector<ComplexMatrix> &wfc_k)
{
    // cout<<file_path<<endl;
    ifstream infile;
    infile.open(file_path);
    string rvalue, ivalue, ik;
    while (infile.peek() != EOF)
    {
        infile >> ik;
        // cout<<"     ik: "<<ik<<endl;
        if (infile.peek() == EOF)
            break;
        // for aims !!!
        for (int iw = 0; iw != NLOCAL; iw++)
            for (int ib = 0; ib != NBANDS; ib++)
                for (int is = 0; is != NSPIN; is++)
                {
                    // cout<<iw<<ib<<is<<ik;
                    infile >> rvalue >> ivalue;
                    // cout<<rvalue<<ivalue<<endl;
                    wfc_k.at(stoi(ik) - 1)(ib, iw) = complex<double>(stod(rvalue), stod(ivalue));
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
void handle_Cs_file(const std::string &file_path)
{
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    natom = stoi(natom_s);
    ncell = stoi(ncell_s);
    // cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl;
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
        if ((*cs_ptr).absmax() >= cs_threshold)
            Cs[ia1][ia2][box] = cs_ptr;

        // cout<<" READ Cs, INDEX:  "<<ia1<<"   "<<ia2<<"   "<<box<<"   "<<(*Cs.at(ia1).at(ia2).at(box))(n_i*n_j-1,n_mu-1)<<endl;
    }
}
void handle_Vq_file(const std::string &file_path, map<Vector3_Order<double>, ComplexMatrix> &Vq_full)
{
    cout << "Begin to read aims vq_real" << endl;
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
        Vector3_Order<double> qvec(kvec_c[stoi(q_num) - 1]);
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

void READ_AIMS_Vq(const std::string &file_path)
{
    const double vq_threshold = 1e-6;
    int vq_save = 0;
    int vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(file_path.c_str());
    vector<string> files;
    map<Vector3_Order<double>, ComplexMatrix> Vq_full;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("coulomb_mat") == 0)
        {
            handle_Vq_file(fm, Vq_full);
        }
    }
    cout << "FINISH coulomb files reading!" << endl;
    for (auto &vf_p : Vq_full)
    {
        auto qvec = vf_p.first;
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
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(J, i_nu), atom_mu_loc2glo(I, i_mu)); ////for aims
                        //(*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                // if (I == J)
                // {
                //     (*vq_ptr).set_as_identity_matrix();
                // }

                if ((*vq_ptr).real().absmax() >= vq_threshold)
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
    cout << "vq threshold: " << vq_threshold << endl;
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

void READ_AIMS_Cs(const std::string &file_path)
{
    cout << "Begin to read Cs" << endl;
    cout << "cs_threshold:  " << cs_threshold << endl;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(file_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("Cs_data") == 0)
            handle_Cs_file(fm);
    }
    closedir(dir);
    dir = NULL;
    cout << "Finish read Cs" << endl;
}

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_nw[ia];
    return iw_lcoal + nb;
}

int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_mu[ia];
    return mu_lcoal + nb;
}

double get_E_min_max(double &Emin, double &Emax)
{
    double E_homo = -1e7;
    double E_lumo = 1e7;
    double midpoint = 1.0 / (NSPIN * n_kpoints);
    double max_lb = ekb[0][0];
    double max_ub = ekb[0][NBANDS - 1];
    cout << "midpoint:  " << midpoint << endl;
    for (int ik = 0; ik != n_kpoints; ik++)
    {
        max_lb = (max_lb > ekb[ik][0]) ? ekb[ik][0] : max_lb;
        max_ub = (max_ub < ekb[ik][NBANDS - 1]) ? ekb[ik][NBANDS - 1] : max_ub;
        int homo_level = 0;
        for (int n = 0; n != NBANDS; n++)
        {
            if (wg(ik, n) >= midpoint)
            {
                homo_level = n;
            }
        }
        cout << "ik  " << ik << "    hommo_level:" << homo_level << endl;
        if (E_homo < ekb[ik][homo_level])
            E_homo = ekb[ik][homo_level];

        if (E_lumo > ekb[ik][homo_level + 1])
            E_lumo = ekb[ik][homo_level + 1];
        // cout<<"    ik:"<<ik<<"     LUMO_level:"<<LUMO_level<<"     E_homo:<<"E_homo<<"      E_lumo:"<<E_lumo<<endl;
    }
    double gap;
    // gap=0.5*(wf.ekb[0][LUMO_level]-wf.ekb[0][LUMO_level-1]);
    Emax = 0.5 * (max_ub - max_lb);
    Emin = 0.5 * (E_lumo - E_homo);
    gap = Emin;
    cout << "E_max: " << Emax << endl;
    cout << "E_min: " << Emin << endl;
    std::cout << "E_homo(ev):" << E_homo * 0.5 * HA2EV << "     E_lumo(ev):" << E_lumo * 0.5 * HA2EV << "     gap(ev): " << gap * HA2EV << endl;
    return gap;
}

