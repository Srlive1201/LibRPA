#ifndef TIMEFREQ_H
#define TIMEFREQ_H

#include <map>
#include <string>
using std::map;
using std::string;

extern const string minimax_grid_path;
extern const string GX_path;

map<double, double> read_local_grid(int grid_N, const string &file_path, const char type, double gap);
#endif
