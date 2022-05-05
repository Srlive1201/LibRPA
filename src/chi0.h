/*!
 @file chi0.h
 @brief Utlities to compute the independent response function
 @warning Not work now
 */
#pragma once
#include <vector>
#include "meanfield.h"
#include "timefreq.h"
#include "ri.h"
#include "vector3_order.h"

using std::vector;

class Chi0
{
    private:
        MeanField mf;
        TFGrids tfgrids;
        double gf_threshold;
        char gf_timefreq;
        char chi_timefreq;
    public:
        Chi0(const MeanField &meanf,
             const TFGrids &tfg,
             const vector<Vector3_Order<int>> &R_list,
             const vector<Vector3_Order<int>> &k_list,
             double threshold = 1.0e-9):
            mf(meanf), tfgrids(tfg), gf_threshold(threshold) {};
        ~Chi0();
        void compute();
};
