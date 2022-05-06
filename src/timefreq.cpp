#include "timefreq.h"
#include "envs.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

using std::map;
using std::pair;
using std::string;
using std::vector;
using std::ifstream;
using std::endl;
using std::cout;

const string minimax_grid_path = string(source_dir) + "/minimax_grid";
const string GX_path = minimax_grid_path + "/GreenX/generate_local_grid.py";

map<double, double> read_local_grid(int grid_N, const string &file_path, const char type, double scale)
{
    ifstream infile;
    map<double, double> grid;
    infile.open(file_path);

    vector<double> tran(grid_N * 2);
    string ss;
    int itran = 0;
    while (infile.peek() != EOF)
    {
        if (infile.peek() == EOF)
            break;
        infile >> ss;
        tran[itran] = stod(ss);
        itran++;
    }
    for (int i = 0; i != grid_N; i++)
    {
        grid.insert(pair<double, double>(tran[i], tran[i + grid_N]));
    }

    infile.close();
    map<double, double> minimax_grid;
    minimax_grid.clear();
    switch (type)
    {
    case 'F':
    {
        if (grid_N <= 20)
        {
            for (auto &i_pair : grid)
                minimax_grid.insert({i_pair.first * scale, i_pair.second * scale * 0.25});
        }
        else
        {
            for (auto &i_pair : grid)
                minimax_grid.insert({i_pair.first * scale, i_pair.second * scale});
        }

        cout << " MINIMAX_GRID_Freq " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }

    case 'T':
    {
        for (auto i_pair : grid)
            minimax_grid.insert({i_pair.first / (scale), i_pair.second / (scale)});

        cout << " MINIMAX_GRID_Tau " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }
    }

    return minimax_grid;
}

vector<double> read_time2freq_trans(const string &file_path, double inverse_scale)
{
    ifstream infile;
    infile.open(file_path);

    vector<double> tran;
    string s;
    int n_freq = 0;
    // stringstream ss;
    // double ss_d;
    while (getline(infile, s))
    {
        if (s[0] == 'F')
        {
            n_freq++;
        }
        else
        {
            stringstream ss(s);
            double ss_d;
            ss >> ss_d;
            tran.push_back(ss_d);
        }
    }
    infile.close();

    double gap;
    double Emin, Emax;
    /* cout << "Cosine_tran_grid" << endl; */
    cout << "read transformation grid" << endl;
    for (int i = 0; i != tran.size(); i++)
    {
        tran[i] /= inverse_scale;
        // cout<<tran[i]<<endl;
    }
    return tran;
}

const string TFGrids::GRID_TYPES_NOTES[TFGrids::N_GRID_TYPES] =
    {
        "Gauss-Legendre grids",
        "Gauss-Chebyshev grids of the first kind",
        "Gauss-Chebyshev grids of the second kind",
        "Minimax grids",
    };

const bool TFGrids::SUPPORT_TIME_GRIDS[TFGrids::N_GRID_TYPES] = 
    { false, false, false, true };

void TFGrids::set_freq()
{
    freq_nodes.resize(n_grids);
    freq_weights.resize(n_grids);
}

void TFGrids::set_time()
{
    time_nodes.resize(n_grids);
    time_weights.resize(n_grids);
    costrans_t2f.create(n_grids, n_grids);
    sintrans_t2f.create(n_grids, n_grids);
    fourier_t2f.create(n_grids, n_grids);
}

void TFGrids::unset()
{
    freq_nodes.clear();
    freq_weights.clear();
    time_nodes.clear();
    time_weights.clear();
    costrans_t2f.create(0, 0);
    sintrans_t2f.create(0, 0);
    fourier_t2f.create(0, 0);
}

TFGrids::TFGrids(unsigned N)
{
    n_grids = N;
    set_freq();
}

void TFGrids::reset(unsigned N)
{
    unset();
    n_grids = N;
    set_freq();
}

TFGrids::~TFGrids()
{
    /* unset(); */
}

void TFGrids::generate_minimax(double emin, double emax)
{
    grid_type = TFGrids::GRID_TYPES::Minimax;
    set_time();
    string tmps;
    // WARN: handle zero division
    double erange = emax / emin;
    tmps = "python " + GX_path + " " + to_string(n_grids) + " " + to_string(erange);
    system(tmps.c_str());

    map<double, double> freq_grid = read_local_grid(n_grids, "local_" + to_string(n_grids) + "_freq_points.dat", 'F', emin);
    int ig = 0;
    for (auto nw: freq_grid)
    {
        freq_nodes[ig] = nw.first;
        freq_weights[ig] = nw.second;
        ig++;
    }
    map<double, double> time_grid = read_local_grid(n_grids, "local_" + to_string(n_grids) + "_time_points.dat", 'T', emin);
    ig = 0;
    for (auto nw: time_grid)
    {
        time_nodes[ig] = nw.first;
        time_weights[ig] = nw.second;
        ig++;
    }
    vector<double> trans;
    // cosine transform
    trans = read_time2freq_trans(to_string(n_grids) + "_time2freq_grid_cos.txt", erange);
    for (int k = 0; k != n_grids; k++)
        for (int j = 0; j != n_grids; j++)
            costrans_t2f(k, j) = trans[ k * n_grids + j];
    trans = read_time2freq_trans(to_string(n_grids) + "_time2freq_grid_sin.txt", erange);
    for (int k = 0; k != n_grids; k++)
        for (int j = 0; j != n_grids; j++)
            sintrans_t2f(k, j) = trans[ k * n_grids + j];
}
