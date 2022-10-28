#include "envs.h"
#include <cctype>
#include <cstring>

//! Get the directory name of absolute file path file_abspath
/*!
 @param[in] file_abspath: absolute path of a file. Support UNIX only.
 */
const char * get_dirname(const char * file_abspath)
{
    char * fn = new char [std::strlen(file_abspath) + 1];
    std::strcpy(fn, file_abspath);
    char *p;
    p = std::strrchr(fn, '/');
    if(p) p[0] = '\0';
    int l = std::strlen(fn);
    char * dn = new char [l + 1];
    std::strcpy(dn, fn);
    return fn;
}

const char * source_dir = get_dirname(__FILE__);
