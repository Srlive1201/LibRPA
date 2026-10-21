#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/timefreq.h"
#include "../io/global_io.h"
#include "../io/stl_io_helper.h"
#include "../mpi/global_mpi.h"
#include "mpi_test_config.h"
#include "testutils.h"

using namespace std;
using namespace librpa_int;

void check_initialize()
{
    cout << "Available time-frequency grids: " << LIBRPA_TFGRID_COUNT << endl;
    // cout << source_dir << endl;
    // cout << minimax_grid_path << endl;
    TFGrids tfg(6);
}

void check_gauss_grids()
{
    TFGrids gl(4);
    gl.generate(LIBRPA_TFGRID_GAUSS_LEGENDRE);
    assert(gl.get_grid_type() == LIBRPA_TFGRID_GAUSS_LEGENDRE);
    assert(!gl.has_time_grids());

    const auto gl_freq = gl.get_freq_nodes();
    const auto gl_weight = gl.get_freq_weights();
    const vector<double> gl_ref_freq = {
        0.037306157410633582,
        0.24627921401387612,
        1.0151079984602933,
        6.7013066301151962,
    };
    const vector<double> gl_ref_weight = {
        0.10042496565842267,
        0.36320093923997282,
        1.4970332756138096,
        18.039340819487776,
    };
    const double tol = 1e-13;
    for (size_t i = 0; i != gl.size(); ++i)
    {
        assert(std::isfinite(gl_freq[i]) && gl_freq[i] > 0.0);
        assert(std::isfinite(gl_weight[i]) && gl_weight[i] > 0.0);
        assert(std::fabs(gl_freq[i] - gl_ref_freq[i]) < tol);
        assert(std::fabs(gl_weight[i] - gl_ref_weight[i]) < tol);
        if (i > 0) assert(gl_freq[i - 1] < gl_freq[i]);
    }

    TFGrids split_gl(4);
    split_gl.generate(LIBRPA_TFGRID_SPLIT_GAUSS_LEGENDRE, -1, -1, 2.0);
    assert(split_gl.get_grid_type() == LIBRPA_TFGRID_SPLIT_GAUSS_LEGENDRE);
    assert(!split_gl.has_time_grids());

    const auto split_freq = split_gl.get_freq_nodes();
    const auto split_weight = split_gl.get_freq_weights();
    const vector<double> split_ref_freq = {
        0.42264973081037427,
        2.2540333075851664,
        4.0,
        17.745966692414836,
    };
    const vector<double> split_ref_weight = {
        1.0,
        0.70564807662546158,
        3.5555555555555558,
        43.738796367818985,
    };
    for (size_t i = 0; i != split_gl.size(); ++i)
    {
        assert(std::isfinite(split_freq[i]) && split_freq[i] > 0.0);
        assert(std::isfinite(split_weight[i]) && split_weight[i] > 0.0);
        assert(std::fabs(split_freq[i] - split_ref_freq[i]) < tol);
        assert(std::fabs(split_weight[i] - split_ref_weight[i]) < tol);
        if (i > 0) assert(split_freq[i - 1] < split_freq[i]);
    }
}

void check_minimax_narrow_transition_interval()
{
    // Actual H2 HF/DZP occupied-to-empty transitions distinguish a usable
    // transform from GreenX success with non-finite fitting coefficients.
    const std::array<double, 9> transitions{
        0.857382174146879, 1.234410186662679, 1.256077372121357,
        1.256077372121389, 2.107728306226503, 2.536985168746743,
        2.593259340292944, 2.593259340293027, 2.877441176675997};
    const auto check_transform = [&](const TFGrids& grid)
    {
        const auto times = grid.get_time_nodes();
        const auto frequencies = grid.get_freq_nodes();
        const auto& transform = grid.get_costrans_t2f();
        for (std::size_t i = 0; i < grid.size(); ++i)
            for (const double delta : transitions)
            {
                double actual = 0.0;
                for (std::size_t j = 0; j < grid.size(); ++j)
                    actual += transform(i, j) * std::exp(-delta * times[j]);
                const double expected =
                    2.0 * delta / (delta * delta + frequencies[i] * frequencies[i]);
                if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-7)
                    throw std::runtime_error("native minimax transition transform is inaccurate");
            }
    };
    TFGrids narrow(24);
    bool rejected = false;
    try
    {
        narrow.generate_minimax(0.4286910870734395, 3.021313235509797);
    }
    catch (const std::exception& error)
    {
        rejected = true;
        cout << "Rejected unusable narrow minimax grid: " << error.what() << endl;
    }
    // An improved GreenX may eventually support this interval without rejecting
    // it; either an accurate transform or rejection before numerical use is valid.
    if (!rejected) check_transform(narrow);
    TFGrids common(24);
    common.generate_minimax(0.3619945912373965, 4.874507376794855);
    check_transform(common);
    narrow.generate_minimax(0.4286910870734395, 3.021313235509797, 1e-12);
    check_transform(narrow);
}

void check_minimax_ng16_diamond_k222()
{
    TFGrids tfg(16);
    // Check minimax grids
    // emin and emax (in Hartree) from diamond PBE k2c
    double emin = 0.173388;
    double emax = 23.7738;
    double residual = tfg.generate_minimax(emin, emax);
    cout << tfg.get_freq_nodes() << endl;
    // TODO(minye): CI failed after this part. I am completely ignorant of why it should.
    cout << tfg.get_freq_weights() << endl;
    cout << tfg.get_time_nodes() << endl;
    cout << tfg.get_time_weights() << endl;
    cout << "Residual: " << residual << endl;
    /* if (tfg.get_grid_type() != LIBRPA_TFGRID_MINIMAX) */
    /*     throw logic_error("internal type should be minimax grid"); */
    /* print_matrix("cos: t2f * f2t, ideally close to identity", tfg.get_costrans_t2f() * tfg.get_costrans_f2t() ); */
    /* print_matrix("sin: t2f * f2t, ideally close to identity", tfg.get_sintrans_t2f() * tfg.get_sintrans_f2t() ); */
    // tfg.write_cos_sin_trans_matrices("test.dat");
}

void check_minimax_ng6_HF_123()
{
    // TODO(minyez): the correct data vary due to the version of GreenX data.
    // The following is correct for the local grid data, but not for the latest GreenX library
    // Need to adapt later.
    TFGrids tfg(6);
    double emin = 0.657768, emax = 30.1366;
    tfg.generate_minimax(emin, emax);
    assert ( tfg.size() == 6 );
    /* printf("%20.15f, %20.15f\n", tfg.get_freq_nodes()[0], tfg.get_freq_weights()[0]); */
    vector<double> freq_node = {0.233556, 0.844872, 2.029850, 4.815547, 12.239097, 36.979336};
    vector<double> freq_weight = {0.489089, 0.798829, 1.724256, 4.282121, 12.092283, 48.530744};
    vector<double> time_node = {0.021614, 0.129568, 0.408121, 1.072294, 2.593619, 6.074920};
    vector<double> time_weight = {0.057049, 0.171324, 0.419863, 0.982601, 2.222569, 5.258784};
    /* cout << tfg.find_freq_weight(tfg.get_freq_nodes()[0]) << endl; */
    assert( 1 == tfg.get_freq_index(tfg.get_freq_nodes()[1]));
    assert( tfg.find_freq_weight(tfg.get_freq_nodes()[0]) == tfg.get_freq_weights()[0]);
    assert( 2 == tfg.get_time_index(tfg.get_time_nodes()[2]));
    auto i_tf = tfg.get_tf_index({tfg.get_time_nodes()[2], tfg.get_freq_nodes()[1]});
    assert(i_tf.first  == tfg.get_time_index(tfg.get_time_nodes()[2]));
    assert(i_tf.second  == tfg.get_time_index(tfg.get_time_nodes()[1]));
    for ( size_t i = 0; i != tfg.size(); i++ )
    {
        assert( fabs(freq_node[i] - tfg.get_freq_nodes()[i]) < 1e-5);
        assert( fabs(freq_weight[i] - tfg.get_freq_weights()[i]) < 1e-5);
        // assert( fabs(time_node[i] - tfg.get_time_nodes()[i]) < 1e-5);
        // assert( fabs(time_weight[i] - tfg.get_time_weights()[i]) < 1e-5);
    }
    // matrix costrans_t2f(6, 6);
    // costrans_t2f(0, 0) = 0.11418;
    // costrans_t2f(0, 1) = 0.34210;
    // costrans_t2f(0, 2) = 0.83715;
    // costrans_t2f(0, 3) = 1.90047;
    // costrans_t2f(0, 4) = 3.66221;
    // costrans_t2f(0, 5) = 1.66271;
    // costrans_t2f(1, 0) = 0.11384;
    // costrans_t2f(1, 1) = 0.34204;
    // costrans_t2f(1, 2) = 0.78564;
    // costrans_t2f(1, 3) = 1.21529;
    // costrans_t2f(1, 4) = -2.37093;
    // costrans_t2f(1, 5) = -2.80551;
    // costrans_t2f(2, 0) = 0.11532;
    // costrans_t2f(2, 1) = 0.32304;
    // costrans_t2f(2, 2) = 0.59313;
    // costrans_t2f(2, 3) = -1.09755;
    // costrans_t2f(2, 4) = -0.39995;
    // costrans_t2f(2, 5) = 2.31135;
    // costrans_t2f(3, 0) = 0.11016;
    // costrans_t2f(3, 1) = 0.29865;
    // costrans_t2f(3, 2) =-0.37431;
    // costrans_t2f(3, 3) =-0.21496;
    // costrans_t2f(3, 4) = 0.52021;
    // costrans_t2f(3, 5) =-1.71560;
    // costrans_t2f(4, 0) = 0.11853;
    // costrans_t2f(4, 1) =-0.05220;
    // costrans_t2f(4, 2) =-0.16181;
    // costrans_t2f(4, 3) = 0.21842;
    // costrans_t2f(4, 4) =-0.36024;
    // costrans_t2f(4, 5) = 1.24194;
    // costrans_t2f(5, 0) = 0.04570;
    // costrans_t2f(5, 1) =-0.08751;
    // costrans_t2f(5, 2) = 0.08781;
    // costrans_t2f(5, 3) =-0.10941;
    // costrans_t2f(5, 4) = 0.19074;
    // costrans_t2f(5, 5) =-0.67673;
    // assert ( is_mat_A_equal_B(6, 6, costrans_t2f.c, tfg.get_costrans_t2f().c, false, false, 1e-5) );
    /* print_matrix("cos: t2f * f2t, ideally close to identity", tfg.get_costrans_t2f() * tfg.get_costrans_f2t() ); */
    /* print_matrix("sin: t2f * f2t, ideally close to identity", tfg.get_sintrans_t2f() * tfg.get_sintrans_f2t() ); */

/* Sine transform matrix */
/*  */
/*   */
/*    0.00072   0.00950   0.08336   0.47439   2.58887   9.87802 */
/*    0.00074   0.04589   0.24892   1.68667   2.90489  -6.70134 */
/*    0.00845   0.06669   0.72122   1.11313  -1.91638   3.95426 */
/*    0.00614   0.24116   0.53245  -0.64767   0.70235  -2.01430 */
/*    0.03966   0.25076  -0.20340   0.15712  -0.22819   0.75918 */
/*    0.06496  -0.00675  -0.01183   0.01955  -0.03648   0.13263 */
}

void check_minimax_ng32_H2O()
{
    // For large minimax grids, high precision is required to reproduce the actual results
    const double emin = 0.302066195818780992e-2 + 0.260399583086844855e0;
    const double emax = 0.322960860707412234e1 + 0.187886173733204913e2;
    // cout << "emin/erange: " << emin << " " << emax/emin << endl;
    double residual;

    TFGrids tfg(32);
    residual = tfg.generate_minimax(emin, emax);
    cout << "Residual: " << residual << endl;
    // tfg.write_cos_sin_trans_matrices("test.dat");
    residual = tfg.generate_minimax(emin, emax, 1e-12);
    cout << "Residual (regulation): " << residual << endl;
    // tfg.write_cos_sin_trans_matrices("test_regulated.dat");
}

int main (int argc, char **argv)
{
    using namespace librpa_int::global;
    int provided;
    MPI_Init_thread(&argc, &argv, LIBRPA_MPI_THREAD_LEVEL, &provided);
    init_global_mpi(MPI_COMM_WORLD);
    init_global_io();

    check_initialize();
    check_gauss_grids();
    check_minimax_narrow_transition_interval();
    check_minimax_ng16_diamond_k222();
    check_minimax_ng6_HF_123();
    check_minimax_ng32_H2O();

    finalize_global_io();
    finalize_global_mpi();
    MPI_Finalize();
    return 0;
}
