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
 */
class TFGrids
{
    public:
        static const int N_GRID_TYPES = 4;
        enum GRID_TYPES { GaussLegendre = 0,
                          GaussChebyshevI = 1,
                          GaussChebyshevII = 2,
                          Minimax = 3, };
        static const string GRID_TYPES_NOTES[N_GRID_TYPES];
        static const bool SUPPORT_TIME_GRIDS[N_GRID_TYPES];
    private:
        //! Internal storage of grid type
        GRID_TYPES grid_type;
        //! whether to use time grids, i.e. space-time method
        bool _use_time_grids;
        size_t n_grids;
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
        //! Parse grid type and the number of grids, do the internal check
        template <typename T>
        void parse_grid_type_n(GRID_TYPES gt, const T &N, bool use_tg);
        //! set the interval of time-freq integration
        /* void set_freq_int_interval(double llim = 0, double ulim = PINF); */
        //! allocate the pointers of array, e.g. nodes and weights
        void set();
        //! delete the pointers
        void unset();
        void generate_minimax(double erange);
    public:
        TFGrids(GRID_TYPES gt, int N, bool use_tg);
        // disable copy at present
        TFGrids(const TFGrids &tfg) {};
        void reset(GRID_TYPES gt, int N, bool use_tg);
        size_t get_n_grids() { return n_grids; }
        vector<double> get_freq_nodes() { return freq_nodes; }
        vector<double> get_freq_weights() { return freq_weights; }
        vector<double> get_time_nodes() { return time_nodes; }
        vector<double> get_time_weights() { return time_weights; }
        //! generate the grid nodes and weights, and transformation matrices if applicable
        void generate();
        void generate(double param1);
        bool use_time_grids() { return _use_time_grids; }
        ~TFGrids();
};

#endif
