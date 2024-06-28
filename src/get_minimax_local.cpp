#include "get_minimax.h"

#include "envs_dir.h"
#include "envs_mpi.h"
#include "parallel_mpi.h"

#include <fstream>
#include <unistd.h>
#include <vector>
#include <string>

static const std::string minimax_grid_path = string(LIBRPA::envs::source_dir) + "/minimax_grid";
static const std::string GX_path = minimax_grid_path + "/GreenX/generate_local_grid.py";

//! read the file containing grids points information
/*!
 @param[in] grid_N: number of grid points
 @param[in] file_path: the file containing the grids
 @param[in] type: the type of grid, either time ('T') or frequency ('F')
 @param[in] scale: scaling parameter
 */
static std::vector<std::pair<double, double>> read_local_grid(int n_grids, const string &file_path, const char type, double scale)
{
    std::ifstream infile;
    std::vector<std::pair<double, double>> grid;
    infile.open(file_path);

    std::vector<double> tran(n_grids * 2);
    std::string ss;
    int itran = 0;
    while (infile.peek() != EOF)
    {
        if (infile.peek() == EOF)
            break;
        infile >> ss;
        tran[itran] = stod(ss);
        itran++;
    }
    for (int i = 0; i != n_grids; i++)
    {
        grid.push_back({tran[i], tran[i + n_grids]});
    }

    infile.close();
    std::vector<std::pair<double, double>> minimax_grid;
    minimax_grid.clear();
    switch (type)
    {
        case 'F':
        {
            if (n_grids <= 20)
            {
                for (auto &i_pair : grid)
                    minimax_grid.push_back({i_pair.first * scale, i_pair.second * scale * 0.25});
            }
            else
            {
                for (auto &i_pair : grid)
                    minimax_grid.push_back({i_pair.first * scale, i_pair.second * scale});
            }
            break;
        }

        case 'T':
        {
            for (auto i_pair : grid)
                minimax_grid.push_back({i_pair.first / (scale), i_pair.second / (scale)});
            break;
        }
    }

    return minimax_grid;
}

//! read the file containing transformation matrix from time to frequency domain
static std::vector<double> read_trans_matrix(const string &file_path, double inverse_scale)
{
    std::ifstream infile;
    infile.open(file_path);

    std::vector<double> tran;
    std::string s;
    // stringstream ss;
    // double ss_d;
    while (getline(infile, s))
    {
        if (s[0] == 'F' || s[0] == 'T')
            continue;
        else
        {
            stringstream ss(s);
            double ss_d;
            ss >> ss_d;
            tran.push_back(ss_d);
        }
    }
    infile.close();

    for (int i = 0; i != tran.size(); i++)
    {
        tran[i] /= inverse_scale;
    }
    return tran;
}

static void call_local_grid_script(int ngrids, double e_min, double e_max)
{
    string tmps;
    double erange = e_max / e_min;
    if(LIBRPA::envs::mpi_comm_global_h.is_root())
    {
        tmps = "python " + GX_path + " " + to_string(ngrids) + " " + to_string(erange);
        system(tmps.c_str());
    }
    LIBRPA::envs::mpi_comm_global_h.barrier();
    sleep(1);
}

static int check_inquiry(int ngrids, double e_min, double e_max)
{
    int ierr = 0;
    return ierr;
}

void get_minimax_grid_frequency(int ngrids, double e_min, double e_max, double *omega_points, double *omega_weights, int &ierr)
{
    ierr = check_inquiry(ngrids, e_min, e_max);
    if (ierr != 0)
        return;
    call_local_grid_script(ngrids, e_min, e_max);
    auto freq_grids = read_local_grid(ngrids, "local_" + to_string(ngrids) + "_freq_points.dat", 'F', e_min);
    for (int i = 0; i < ngrids; i++)
    {
        omega_points[i] = freq_grids[i].first;
        omega_weights[i] = freq_grids[i].second;
    }
}

void get_minimax_grid(int ngrids, double e_min, double e_max,
                      double *tau_points, double *tau_weights,
                      double *omega_points, double *omega_weights,
                      double *cosft_wt, double *cosft_tw, double *sinft_wt,
                      double max_errors[3], double &cosft_duality_error,
                      int &ierr)
{
    ierr = check_inquiry(ngrids, e_min, e_max);
    if (ierr != 0)
        return;
    call_local_grid_script(ngrids, e_min, e_max);
    auto freq_grids = read_local_grid(ngrids, "local_" + to_string(ngrids) + "_freq_points.dat", 'F', e_min);
    auto time_grids = read_local_grid(ngrids, "local_" + to_string(ngrids) + "_time_points.dat", 'T', e_min);
    for (int i = 0; i < ngrids; i++)
    {
        omega_points[i] = freq_grids[i].first;
        omega_weights[i] = freq_grids[i].second;
        tau_points[i] = time_grids[i].first;
        tau_weights[i] = time_grids[i].second;
    }

    vector<double> trans;
    trans = read_trans_matrix(to_string(ngrids) + "_time2freq_grid_cos.txt", e_min);
    for (int i = 0; i < ngrids * ngrids; i++)
        cosft_wt[i] = trans[i];
    trans = read_trans_matrix(to_string(ngrids) + "_freq2time_grid_cos.txt", 1/e_min);
    for (int i = 0; i < ngrids * ngrids; i++)
        cosft_tw[i] = trans[i];
    trans = read_trans_matrix(to_string(ngrids) + "_time2freq_grid_sin.txt", e_min);
    for (int i = 0; i < ngrids * ngrids; i++)
        sinft_wt[i] = trans[i];

    // TODO(minyez): read the error from the script log

    // since the matrix is very small (rows/cols<100), so we do manual matmul here
    cosft_duality_error = -1.0;
    for (int i = 0; i < ngrids; i++)
        for (int j = 0; j < ngrids; j++)
        {
            double error_ij = i == j ? -1.0: 0.0;
            for (int k = 0; k < ngrids; k++)
                error_ij += cosft_wt[i*ngrids + k] * cosft_tw[k*ngrids + j];
            cosft_duality_error = std::max(cosft_duality_error, std::abs(error_ij));
        }
}
