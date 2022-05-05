#include "input.h"
/* #include "atoms.h" */
#include <fstream>
#include <iostream>
#include "cal_periodic_chi0.h"
#include "meanfield.h"
#include "constants.h"
using namespace std;
double cs_threshold = 1E-6;
double vq_threshold = 1e-6;

int kv_nmp[3] = {1, 1, 1};
Vector3<double> *kvec_c;
Matrix3 latvec;
Matrix3 G;

void READ_AIMS_STRU(const std::string &file_path)
{
    int n_kpoints = meanfield.get_n_kpoints();
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
        cout << "kvec [" << i << "]: " << kvec_c[i] << endl;
    }
}
