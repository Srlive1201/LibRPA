#include "timefreq.h"
#include "envs.h"
#include "mathtools.h"
#include "parallel_mpi.h"
#include "get_minimax.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <unistd.h>

using std::pair;
using std::string;
using std::ifstream;
using std::endl;
using std::cout;

const string TFGrids::GRID_TYPES_NOTES[TFGrids::GRID_TYPES::COUNT] =
    {
        "Gauss-Legendre grids",
        "Gauss-Chebyshev grids of the first kind",
        "Gauss-Chebyshev grids of the second kind",
        "Minimax time-frequency grids",
        "Even-spaced frequency grids",
        "Even-spaced time-frequency grids (debug use)",
    };

const bool TFGrids::SUPPORT_TIME_GRIDS[TFGrids::GRID_TYPES::COUNT] =
    { false, false, false, true, false, true };

TFGrids::GRID_TYPES TFGrids::get_grid_type(const string& grid_str)
{
    if (grid_str == "GL")
        return TFGrids::GRID_TYPES::GaussLegendre;
    if (grid_str == "GC-I")
        return TFGrids::GRID_TYPES::GaussChebyshevI;
    if (grid_str == "GL-II")
        return TFGrids::GRID_TYPES::GaussChebyshevII;
    if (grid_str == "minimax")
        return TFGrids::GRID_TYPES::Minimax;
    if (grid_str == "evenspaced")
        return TFGrids::GRID_TYPES::EvenSpaced;
    if (grid_str == "evenspaced_tf")
        return TFGrids::GRID_TYPES::EvenSpaced_TF;
    throw std::invalid_argument("Unknown grid string: " + grid_str);
}

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
    costrans_f2t.create(n_grids, n_grids);
    sintrans_f2t.create(n_grids, n_grids);
    fourier_t2f.create(n_grids, n_grids);
}

void TFGrids::show()
{
    cout << "Grid type: " << TFGrids::GRID_TYPES_NOTES[grid_type] << endl;
    cout << "Grid size: " << n_grids << endl;
    cout << "Frequency node & weight: " << endl;
    for ( int i = 0; i != n_grids; i++ )
        printf("%2d %20.12f %20.12f\n", i, freq_nodes[i], freq_weights[i]);
    if (has_time_grids())
    {
        cout << "Time node & weight: " << endl;
        for ( int i = 0; i != n_grids; i++ )
            printf("%2d %20.12f %20.12f\n", i, time_nodes[i], time_weights[i]);
        // cout << "t->f transform: " << endl;
        // if (costrans_t2f.size)
        // {
        //     print_matrix("Cosine transform matrix", costrans_t2f);
        // }
        // if (sintrans_t2f.size)
        // {
        //     print_matrix("Sine transform matrix", sintrans_t2f);
        // }
    }
    printf("\n");
}

void TFGrids::unset()
{
    freq_nodes.clear();
    freq_weights.clear();
    time_nodes.clear();
    time_weights.clear();
    costrans_t2f.create(0, 0);
    sintrans_t2f.create(0, 0);
    costrans_f2t.create(0, 0);
    sintrans_f2t.create(0, 0);
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

void TFGrids::generate_evenspaced(double emin, double interval)
{
    if ( emin <= 0 )
        throw invalid_argument("emin must be positive");
    if ( interval < 0 )
        throw invalid_argument("emin must be non-negative");
    double weight = 1.0 / n_grids;
    for ( int i = 0; i != n_grids; i++)
    {
        freq_nodes[i] = emin + interval * i;
        freq_weights[i] = weight;
    }
}

void TFGrids::generate_evenspaced_tf(double emin, double eintv, double tmin, double tintv)
{
    generate_evenspaced(emin, eintv);
    set_time();
    if ( tmin <= 0 )
        throw invalid_argument("tmin must be positive");
    if ( tintv < 0 )
        throw invalid_argument("tintv must be non-negative");
    double weight = 1.0 / n_grids;
    for ( int i = 0; i != n_grids; i++)
    {
        time_nodes[i] = tmin + tintv * i;
        time_weights[i] = weight;
        // WARN: fake transform matrices
        costrans_t2f(i, i) = weight;
        sintrans_t2f(i, i) = weight;
        costrans_f2t(i, i) = 1/weight;
        sintrans_f2t(i, i) = 1/weight;
    }
}

void TFGrids::generate_minimax(double emin, double emax)
{
    grid_type = TFGrids::GRID_TYPES::Minimax;
    set_time();

    double * omega_points = new double [n_grids];
    double * tau_points = new double [n_grids];
    double * omega_weights = new double [n_grids];
    double * tau_weights = new double [n_grids];
    double max_errors[3];
    double cosft_duality_error;
    int ierr;

    get_minimax_grid(n_grids, emin, emax, tau_points, tau_weights, omega_points, omega_weights,
                     costrans_t2f.c, costrans_f2t.c, sintrans_t2f.c, max_errors, cosft_duality_error, ierr);

    if (ierr != 0)
        throw invalid_argument(string("minimax grids failed, return code: ") + to_string(ierr));
    printf("Cosine transform duality error: %20.12f\n", cosft_duality_error);

    for (int ig = 0; ig != n_grids; ig++)
    {
        freq_nodes[ig] = omega_points[ig];
        freq_weights[ig] = omega_weights[ig];
        time_nodes[ig] = tau_points[ig];
        time_weights[ig] = tau_weights[ig];
    }

    delete [] omega_points;
    delete [] omega_weights;
    delete [] tau_points;
    delete [] tau_weights;

    // zmy debug
    // cout << "Cos transform time -> freq (freq each row)" << endl;
    // cout << costrans_t2f << endl;
    // cout << "Cos transform freq -> time (time each row)" << endl;
    // cout << costrans_f2t << endl;
    // cout << "Cos transform Delta" << endl;
    // cout << costrans_t2f * costrans_f2t << endl;
    // cout << "Sin transform time -> freq (freq each row)" << endl;
    // cout << sintrans_t2f << endl;
    // cout << "Sin transform freq -> time (time each row)" << endl;
    // cout << sintrans_f2t << endl;
    // cout << "Sin transform Delta" << endl;
    // cout << sintrans_t2f * sintrans_f2t << endl;
}

void TFGrids::generate_GaussChebyshevI()
{
    grid_type = TFGrids::GRID_TYPES::GaussChebyshevI;
    double nodes[n_grids], weights[n_grids];
    GaussChebyshevI_unit(n_grids, nodes, weights);
    // transform from [-1,1] to [0, infinity]
    transform_GaussQuad_unit2x0inf(0.0, n_grids, nodes, weights);
    for ( int i = 0; i != n_grids; i++ )
    {
        freq_nodes[i] = nodes[i];
        freq_weights[i] = weights[i];
    }
}

void TFGrids::generate_GaussChebyshevII()
{
    grid_type = TFGrids::GRID_TYPES::GaussChebyshevII;
    double nodes[n_grids], weights[n_grids];
    GaussChebyshevII_unit(n_grids, nodes, weights);
    // transform from [-1,1] to [0, infinity]
    transform_GaussQuad_unit2x0inf(0.0, n_grids, nodes, weights);
    for ( int i = 0; i != n_grids; i++ )
    {
        freq_nodes[i] = nodes[i];
        freq_weights[i] = weights[i];
    }
}

void TFGrids::generate_GaussLegendre()
{
    grid_type = TFGrids::GRID_TYPES::GaussLegendre;
    double nodes[n_grids], weights[n_grids];
    GaussLegendre_unit(n_grids, nodes, weights);
    // transform from [-1,1] to [0, infinity]
    transform_GaussQuad_unit2x0inf(0.0, n_grids, nodes, weights);
    for ( int i = 0; i != n_grids; i++ )
    {
        freq_nodes[i] = nodes[i];
        freq_weights[i] = weights[i];
    }
}

int TFGrids::get_time_index(const double &time) const
{
    if (!has_time_grids())
        throw logic_error("time grids not available");
    auto itr = std::find(time_nodes.cbegin(), time_nodes.cend(), time);
    if ( itr == time_nodes.cend() )
        throw invalid_argument("time not found");
    int i = std::distance(time_nodes.cbegin(), itr);
    return i;
}

int TFGrids::get_freq_index(const double &freq) const
{
    auto itr = std::find(freq_nodes.cbegin(), freq_nodes.cend(), freq);
    if ( itr == freq_nodes.cend() )
        throw invalid_argument("frequency not found");
    int i = std::distance(freq_nodes.cbegin(), itr);
    return i;
}

const pair<int, int> TFGrids::get_tf_index(const pair<double, double> &tf) const
{
    auto itime = get_time_index(tf.first);
    auto ifreq = get_freq_index(tf.second);
    return pair<int, int>{itime, ifreq};
}

double TFGrids::find_freq_weight(const double & freq) const
{
    return freq_weights[get_freq_index(freq)];
}
