#include "fs.h"

#include <filesystem>
#include <stdexcept>
#include <system_error>

#include "../utils/error.h"

namespace librpa_int
{

std::string path_as_directory(const std::string &path)
{
    if (path.find(":") != std::string::npos)
    {
        throw std::runtime_error("dirpath contains invalid character (:) for POSIX path");
    }

    if (path.back() != '/')
    {
        return path + '/';
    }

    return path;
}

void create_directories(const char *dname, int root_process)
{
    if (std::filesystem::is_directory(dname)) return;
    if (root_process == 0)
    {
        std::error_code ec;
        std::filesystem::create_directories(dname, ec);
        if (!std::filesystem::is_directory(dname))
        {
            throw LIBRPA_RUNTIME_ERROR(std::string("Failed to create directories ") + dname);
        }
    }
}

}
