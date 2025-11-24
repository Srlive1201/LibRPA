#pragma once
#include <fstream>

namespace LIBRPA
{

namespace envs
{

extern std::ofstream ofs;

//! File output stream handler of each process
extern std::ofstream ofs_myid;

//! Control whether to redirect stdout (cout, fmt print) to file
extern bool redirect_stdout;

//! File stream used by fprintf when stdout is redirected
extern FILE *pfile_redirect;

//! Initialize the IO environment of LibRPA
/*!
 * @param  [in]  redirect_stdout_in    flag to control to redirect stdout to file, useful when separating LibRPA output from the main program
 * @pram   [in]  filename              name of file to redirect stdout
 */
void initialize_io(bool redirect_stdout_in = false, const char *filename = "LibRPA_output.txt");

//! Check the IO environment of LibRPA is correctly initialized
bool is_io_initialized();


//! Finalize the IO environment of LibRPA
void finalize_io();

} /* end of namespace envs */


} /* end of namespace LIBRPA */
