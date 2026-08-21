#pragma once
/**
 * @file librpa_file_reader.hpp
 * @brief File input interface for LibRPA datasets.
 */

// Public API headers
#include "librpa.hpp"

// Internal data types exposed by this C++ interface
#include "../src/api/dataset.h"
#include "../src/core/coulmat.h"
#include "../src/math/lapack_connector.h"
#include "../src/math/scalapack_connector.h"

// Standard headers
#include <memory>
#include <string>

namespace librpa
{

/** @brief Options controlling dataset input from LibRPA files. */
struct FileReaderOptions
{
    std::string input_dir;               ///< Directory containing LibRPA input files.
    double cs_threshold = 1.0e-12;       ///< Screening threshold for RI coefficients.
    double coulomb_threshold = 1.0e-12;  ///< Screening threshold for Coulomb matrices.
    bool read_ri = true;                 ///< Whether to read RI and Coulomb data.
    bool read_band_data = true;          ///< Whether to read a separate band-path mean field.
};

/**
 * @brief Read a LibRPA dataset from files.
 *
 * @param[in] comm MPI communicator used to initialize the dataset.
 * @param[in] opts File reader options.
 * @return Shared pointer to the populated LibRPA dataset.
 */
std::shared_ptr<librpa_int::Dataset> read_dataset_from_files(
    MPI_Comm comm, const FileReaderOptions &opts);

} // namespace librpa
