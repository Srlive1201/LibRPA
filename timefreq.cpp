#include "timefreq.h"
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

// TODO can minimax_grid_path be a dynamic variable determined at runtime?
const string minimax_grid_path = "/home/minyez/projects/LibRPA/minimax_grid";
const string GX_path = minimax_grid_path + "/GreenX/generate_local_grid.py";

map<double, double> read_local_grid(int grid_N, const string &file_path, const char type, double gap)
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
                minimax_grid.insert({i_pair.first * gap, i_pair.second * gap * 0.25});
        }
        else
        {
            for (auto &i_pair : grid)
                minimax_grid.insert({i_pair.first * gap, i_pair.second * gap});
        }

        cout << " MINIMAX_GRID_Freq " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }

    case 'T':
    {
        for (auto i_pair : grid)
            minimax_grid.insert({i_pair.first / (gap), i_pair.second / (gap)});

        cout << " MINIMAX_GRID_Tau " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }
    }

    return minimax_grid;
}
