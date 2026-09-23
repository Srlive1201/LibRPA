#ifndef READER_BASIS_H
#define READER_BASIS_H

#include <string>

#include "reader_context.h"
namespace librpa::reader
{

void reader_basis(ReaderContext &ctx, const std::string &file_path);
void reader_basis_wfc(ReaderContext &ctx, const std::string &file_path);
void reader_basis_aux(ReaderContext &ctx, const std::string &file_path);
void reader_basis_aux_shrink(ReaderContext &ctx, const std::string &file_path);

}
#endif // READER_BASIS_H
