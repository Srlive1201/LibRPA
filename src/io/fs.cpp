#include "fs.h"
#include <filesystem>
#include <stdexcept>
#include <system_error>

namespace librpa_int
{

void create_directories(const char *dname, int root_process)
{
    if (std::filesystem::is_directory(dname)) return;
    if (root_process == 0)
    {
        std::error_code ec;
        std::filesystem::create_directories(dname, ec);
        if (std::filesystem::is_directory(dname))
        {
            throw std::runtime_error(std::string("Failed to create directories ") + dname);
        }
    }
}

}
