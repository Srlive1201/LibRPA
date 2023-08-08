#include "read_aims.h"
// #include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <algorithm>
#include <parallel_mpi.h>
#include <malloc.h>
#include "atoms.h"
#include "ri.h"
#include "input.h"

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
    }}

size_t READ_AIMS_Cs(const string &dir_path, double threshold)
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
            cs_discard += handle_Cs_file(fm, threshold);
    }
    closedir(dir);
    dir = NULL;
    
    atom_mu_part_range.resize(atom_mu.size());
    atom_mu_part_range[0]=0;
    for(int I=1;I!=atom_mu.size();I++)
        atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // for(int i=0;i!=atom_mu_part_range.size();i++)
    //     cout<<" atom_mu_part_range ,i: "<<i<<"    "<<atom_mu_part_range[i]<<endl;

    // cout << "Finish read Cs" << endl;
    return cs_discard;
}

size_t handle_Cs_file(const string &file_path, double threshold)
{
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

        // if (box == Vector3_Order<int>({0, 0, 1}))continue;
        if ((*cs_ptr).absmax() >= threshold)
            Cs[ia1][ia2][box] = cs_ptr;
        else
            cs_discard++;
        // cout<<" READ Cs, INDEX:  "<<ia1<<"   "<<ia2<<"   "<<box<<"   "<<(*Cs.at(ia1).at(ia2).at(box))(n_i*n_j-1,n_mu-1)<<endl;
    }
    return cs_discard;
}



size_t READ_AIMS_d_Cs(const string &dir_path, double threshold)
{
    size_t cs_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find("d_Cs_data") == 0)
            cs_discard += handle_d_Cs_file(fm, threshold);
    }
    closedir(dir);
    dir = NULL;
    
    atom_mu_part_range.resize(atom_mu.size());
    atom_mu_part_range[0]=0;
    for(int I=1;I!=atom_mu.size();I++)
        atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    return cs_discard;
}

size_t handle_d_Cs_file(const string &file_name, double threshold)
{
    size_t d_cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s;
    string d_Cs_ele_x, d_Cs_ele_y, d_Cs_ele_z;
    ifstream infile;
    infile.open(file_name);
    infile >> natom_s >> ncell_s;
    natom = stoi(natom_s);
    ncell = stoi(ncell_s);
    while (infile.peek() != EOF)   // EOF (end of file reached) 
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF)
            break;
        infile >> j_s >> mu_s;
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

        shared_ptr<matrix> cs_ptr_x = make_shared<matrix>();
        shared_ptr<matrix> cs_ptr_y = make_shared<matrix>();
        shared_ptr<matrix> cs_ptr_z = make_shared<matrix>();

        cs_ptr_x->create(n_i * n_j, n_mu);
        cs_ptr_y->create(n_i * n_j, n_mu);
        cs_ptr_z->create(n_i * n_j, n_mu);

        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> d_Cs_ele_x >> d_Cs_ele_y >> d_Cs_ele_z;

                    (*cs_ptr_x)(i * n_j + j, mu) = stod(d_Cs_ele_x);
                    (*cs_ptr_y)(i * n_j + j, mu) = stod(d_Cs_ele_y);
                    (*cs_ptr_z)(i * n_j + j, mu) = stod(d_Cs_ele_z);
                }

        if ((*cs_ptr_x).absmax() >= threshold)
            d_Cs[ia1][ia2][1][box] = cs_ptr_x;
        else
            d_cs_discard++;


        if ((*cs_ptr_y).absmax() >= threshold)
            d_Cs[ia1][ia2][2][box] = cs_ptr_y;
        else
            d_cs_discard++;

        if ((*cs_ptr_z).absmax() >= threshold)
            d_Cs[ia1][ia2][3][box] = cs_ptr_z;
        else
            d_cs_discard++;

    }
    return d_cs_discard;
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
            handle_Vq_full_file(fm, threshold, Vq_full);
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

void handle_Vq_full_file(const string &file_path, double threshold, map<Vector3_Order<double>, ComplexMatrix> &Vq_full)
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

    while (infile.peek() != EOF)
    {
        infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col;
        if (infile.peek() == EOF)
            break;
        //cout << "vq range: " << begin_row << " ~ " << end_row << "  ,   " << begin_col << " ~ " << end_col << endl;
        infile >> q_num >> q_weight;
        int mu = stoi(nbasbas);
        int nu = stoi(nbasbas);
        int brow = stoi(begin_row) - 1;
        int erow = stoi(end_row) - 1;
        int bcol = stoi(begin_col) - 1;
        int ecol = stoi(end_col) - 1;
        int iq = stoi(q_num) - 1;
        if((erow-brow<=0) || (ecol-bcol<=0))
            return;
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
}


size_t READ_Wc_times_polar(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_mat_t & Wc_times_polar)
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
            handle_Vq_full_file(fm, threshold, Vq_full);
        }
    }
    for (auto &vf_p : Vq_full)
    {
        auto qvec = vf_p.first;
        for (int I = 0; I != atom_mu.size(); I++)
            for (int J = 0; J != atom_mu.size(); J++)
            {
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);
                for (int i_mu = 0; i_mu != atom_mu[I]; i_mu++)
                {

                    for (int i_nu = 0; i_nu != atom_mu[J]; i_nu++)
                    {
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                if ((*vq_ptr).real().absmax() >= threshold)
                {
                    Wc_times_polar[I][J][qvec] = vq_ptr;
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
    return vq_discard;
}



size_t READ_d_Vq_Full(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_d_mat_t &d_coulomb_mat)
{
    size_t vq_save = 0;
    size_t vq_discard = 0;
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(dir_path.c_str());
    vector<string> files;
    map<int,map<int,map<Vector3_Order<double>, ComplexMatrix>>> d_Vq_full;

    while ((ptr = readdir(dir)) != NULL)
    {
        string fm(ptr->d_name);
        if (fm.find(vq_fprefix) == 0)
        {
            handle_d_Vq_full_file(fm, threshold, d_Vq_full);
        }
    }
   
    //for (auto &vf_p : d_Vq_full)
    for (auto &dvq_coord_atom : d_Vq_full)
    {
        auto i_coord = dvq_coord_atom.first;
    for (auto &dvq_atom : dvq_coord_atom.second)
    {	    
        auto i_atom = dvq_atom.first;
    for (auto &vf_p : dvq_atom.second)
    {	    
        auto qvec = vf_p.first;
        for (int I = 0; I != atom_mu.size(); I++)
            for (int J = 0; J != atom_mu.size(); J++)
            {
                if (I > J)
                    continue;
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);

                for (int i_mu = 0; i_mu != atom_mu[I]; i_mu++)
                {

                    for (int i_nu = 0; i_nu != atom_mu[J]; i_nu++)
                    {
//			    cout << vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu))<< "d_coul"<< endl; 
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                if ((*vq_ptr).real().absmax() >= threshold)
                {
                    d_coulomb_mat[i_coord][i_atom][I][J][qvec] = vq_ptr;
		    for(int ir=0; ir<atom_mu[I]; ++ir)
			    for(int ic=0; ic<atom_mu[J]; ++ic)
//		   cout << (*d_coulomb_mat[i_coord][i_atom][I][J][qvec])(ir,ic) << "d_coulomb_mat" <<endl; 
                    vq_save++;
                }
                else
                {
                    vq_discard++;
                }
            }
    }
    }
    }
    closedir(dir);
    dir = NULL;
    return vq_discard;
}


void handle_d_Vq_full_file(const string &file_name, double threshold, map<int, map<int,map<Vector3_Order<double>, ComplexMatrix>>> &d_Vq_full)
{
    ifstream infile;
    infile.open(file_name);
    string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3;
    string vqx_r, vqx_i, vqy_r, vqy_i, vqz_r,vqz_i, q_num, q_weight, natoms;
    infile >> n_irk_points;
    while (infile.peek() != EOF)
    {
        infile >> nbasbas >> begin_row >> end_row >> begin_col >> end_col >> natoms;
        if (infile.peek() == EOF)
            break;
        infile >> q_num >> q_weight;
        int mu = stoi(nbasbas);
        int nu = stoi(nbasbas);
        int brow = stoi(begin_row) - 1;
        int erow = stoi(end_row) - 1;
        int bcol = stoi(begin_col) - 1;
        int ecol = stoi(end_col) - 1;
        int iq = stoi(q_num) - 1;
	int my_n_atoms=stoi(natoms);
        Vector3_Order<double> qvec(kvec_c[iq]);
        if (irk_weight.count(qvec) == 0)
            irk_weight.insert(pair<Vector3_Order<double>, double>(qvec, stod(q_weight)));
	for (int i_atom =1; i_atom<=my_n_atoms; i_atom++)
	{
        for (int i=1; i<=3; i++)
        {			
        if (!d_Vq_full[i][i_atom].count(qvec))
        {
            d_Vq_full[i][i_atom][qvec].create(mu, nu);
        }
	}
//        if (!d_Vq_full[2][i_atom].count(qvec))
//        {
//            d_Vq_full[2][i_atom][qvec].create(mu, nu);
//        }
//        if (!d_Vq_full[3][i_atom].count(qvec))
//        {
//            d_Vq_full[3][i_atom][qvec].create(mu, nu);
//        }
	}
	for (int i_atom =1; i_atom<=my_n_atoms; i_atom++)
        for (int i_mu = brow; i_mu <= erow; i_mu++)
	{
            for (int i_nu = bcol; i_nu <= ecol; i_nu++)
            {
                infile >> vqx_r >> vqx_i >> vqy_r >> vqy_i >> vqz_r >> vqz_i;
	//	cout << vqx_r << vqx_i << " vqx_r;, vqx_i" << endl;
                d_Vq_full[1][i_atom][qvec](i_mu, i_nu) = complex<double>(stod(vqx_r), stod(vqx_i)); 
                d_Vq_full[2][i_atom][qvec](i_mu, i_nu) = complex<double>(stod(vqy_r), stod(vqy_i)); 
                d_Vq_full[3][i_atom][qvec](i_mu, i_nu) = complex<double>(stod(vqz_r), stod(vqz_i)); 
	//   cout << d_Vq_full[3][i_atom][qvec](i_mu, i_nu) << endl;
            }
	}
    }
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

void handle_Vq_row_file(const string &file_path, double threshold, atpair_k_cplx_mat_t &coulomb, const vector<atpair_t> &local_atpair)
{
    // cout << "Begin to read aims vq_real from " << file_path << endl;
    ifstream infile;
    infile.open(file_path);
    string nbasbas, begin_row, end_row, begin_col, end_col, q1, q2, q3, vq_r, vq_i, q_num, q_weight;
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
    for(auto &Ip:Cs)
        if(!loc_atp_index.count(Ip.first))
        {
            Cs.erase(Ip.first);
            
        }
    malloc_trim(0);
    printf("  |process %d, size of Cs after erase: %d,  max_size: %d\n",para_mpi.get_myid(),Cs.size(),Cs.max_size());

}
// TODO: implement the wrapper of all input readers
void read_aims(MeanField &mf)
{
}
