#include "read_aims.h"
#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
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
    }
}

// TODO: implement the wrapper of all input readers
void read_aims(MeanField &mf)
{
}
