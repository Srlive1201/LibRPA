#include "utils_mem.h"

#include <cstdlib>
#if defined(__linux__)|| defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <malloc.h>
#elif defined(__APPLE__) && defined(__MACH__)
#include <malloc/malloc.h>  // for malloc_zone_pressure_relief and malloc_default_zone
#endif

namespace LIBRPA
{

namespace utils
{

void display_free_mem()
{
#if defined(__linux__)
    std::system("free -m");
#elif defined(__APPLE__) && defined(__MACH__)
    std::system("vm_stat | awk 'NR==1{page_size=$8} NR>1{gsub(\".\", \"\"); free+=$3} END {print \"Free memory (in MB): \" page_size*4096/1024/1024}'");
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
    std::system("sysctl hw.physmem hw.usermem");
#endif
}

void release_free_mem()
{
#if defined(__linux__)|| defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
    malloc_trim(0);
#elif defined(__APPLE__) && defined(__MACH__)
    malloc_zone_pressure_relief(malloc_default_zone(), 0);
#endif
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */

