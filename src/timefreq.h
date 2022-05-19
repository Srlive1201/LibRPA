/*!
 @file: timefreq.h
 @brief Utilities related to quadrature on time/frequency domain
 */
#ifndef TIMEFREQ_H
#define TIMEFREQ_H

#include <map>
#include <vector>
#include <string>
#include "constants.h"
#include "matrix.h"
#include "complexmatrix.h"
using std::map;
using std::vector;
using std::string;

extern const string minimax_grid_path;
extern const string GX_path;

//! read the file containing grids points information
/*!
 @param[in] grid_N: number of grid points
 @param[in] file_path: the file containing the grids
 @param[in] type: the type of grid, either time ('T') or frequency ('F')
 @param[in] scale: scaling parameter. Usually the spectral width of energy transition
 */
map<double, double> read_local_grid(int grid_N, const string &file_path, const char type, double scale);
//! read the file containing transformation matrix from time to frequency domain
vector<double> read_time2freq_trans(const string &file_path, double inverse_scale);

//! Object to handle time/frequency grids for quadrature
/*!
 Not necessay have a time grids, unless the space-time minimax grid is used.
 @note Only pure imaginary grid method is implemented.
 @note The time grids are always generated when available, since they are cheap
 */
class TFGrids
{
    public:
        enum GRID_TYPES { GaussLegendre,
                          GaussChebyshevI,
                          GaussChebyshevII,
                          Minimax,
                          EvenSpaced, EvenSpaced_TF,
                          COUNT, // NOTE: always the last
        };
        static const string GRID_TYPES_NOTES[GRID_TYPES::COUNT];
        static const bool SUPPORT_TIME_GRIDS[GRID_TYPES::COUNT];
    private:
        //! Internal storage of grid type
        GRID_TYPES grid_type;
        //! whether to use time grids, i.e. space-time method
        bool _has_time_grids;
        unsigned n_grids;
        vector<double> freq_nodes;
        vector<double> freq_weights;
        vector<double> time_nodes;
        vector<double> time_weights;
        /*! Cosine transformation matrix from (imaginary) time to frequency, i.e. gamma in Eq.4 of LiuP16
         The row-column convention also follows the equation, i.e.
         Each row corresponds to a certain frequency k index.
         */
        matrix costrans_t2f;
        //! Sine transformation matrix from (imaginary) time to frequency, i.e. lambda in Eq.7 of LiuP16
        matrix sintrans_t2f;
        //! General Fourier transformation matrix.
        ComplexMatrix fourier_t2f; // FIXME: this is not implemented, please do not use
        //! allocate the pointers of array, e.g. nodes and weights
        void set_freq();
        void set_time();
        //! delete the pointers
        void unset();
    public:
        TFGrids::GRID_TYPES get_grid_type() const { return grid_type; }
        TFGrids(unsigned N);
        // disable copy at present
        TFGrids(const TFGrids &tfg) {};
        void reset(unsigned N);
        void show();
        //! get the number of grid points
        size_t get_n_grids() const { return n_grids; }
        //! alias to get_n_grids
        size_t size() const { return n_grids; }
        const vector<double> get_freq_nodes() const { return freq_nodes; }
        const vector<double> get_freq_weights() const { return freq_weights; }
        const vector<double> get_time_nodes() const { return time_nodes; }
        const vector<double> get_time_weights() const { return time_weights; }
        // NOTE: use reference for return matrix?
        const matrix get_costrans_t2f() const { return costrans_t2f; }
        const matrix get_sintrans_t2f() const { return sintrans_t2f; }
        //! obtain the integral weight from the frequency value
        // NOTE:(ZMY) attempt to use a map<double, double> to store,
        //      but will lead to a segfault in chi0tauR calculation, not knowing why
        double find_freq_weight(const double &freq);
        //! Generate the even-spaced frequency grid
        void generate_evenspaced(double emin, double interval);
        //! Generate the even-spaced time-frequency grid. @note Currently only for debug use
        void generate_evenspaced_tf(double emin, double eintv, double tmin, double tintv);
        //! Generate the minimax time-frequency grid
        void generate_minimax(double emin, double emax);
        //! Generate Gauss-Chebyshev quadrature of first kind on [0, infty)
        void generate_GaussChebyshevI();
        //! Generate Gauss-Chebyshev quadrature of second kind on [0, infty)
        void generate_GaussChebyshevII();
        void generate_GaussLegendre();
        bool has_time_grids() const { return time_nodes.size() > 0; }
        ~TFGrids();
};

#endif
