#pragma once
#include "../../src/math/matrix_m.h"
#include <map>
#include <string>
namespace librpa::reader {
using VxcMatrices = std::map<int,std::map<int,librpa_int::Matz>>;
// Parse producer Vxc matrices and convert energies to Hartree. The numerical
// task projects AO matrices and expands paired IBZ source states.
VxcMatrices read_qsgw_vxc(const std::string& input_directory,const std::string& prefix,
                          int nspins,int nkpoints,int nspinor,bool aims_input,bool band_path);
}
