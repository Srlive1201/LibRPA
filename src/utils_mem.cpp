#include "utils_mem.h"

#include <cstdlib>
// For memory cleanup
#if defined(__linux__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <malloc.h>
#elif defined(__APPLE__) && defined(__MACH__)
#include <malloc/malloc.h>  // for malloc_zone_pressure_relief and malloc_default_zone
#endif

// For memory query
#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#include <sys/types.h>
#include <sys/sysctl.h>
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
#if defined(__linux__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
    malloc_trim(0);
#elif defined(__APPLE__) && defined(__MACH__)
    malloc_zone_pressure_relief(malloc_default_zone(), 0);
#endif
}

int get_node_total_mem(double &total_mem)
{
    int retcode = 1;
    // value in KB unit
    long value = 0L;

#if defined(__linux__)
    struct sysinfo mem_info;
    retcode = sysinfo(&mem_info);
    value = (mem_info.totalram/1024)*mem_info.mem_unit;
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
    int mib[2];
    size_t len;
    mib[0] = CTL_HW;
    mib[1] = HW_PHYSMEM;
    len = sizeof(value);
    retcode = sysctl(mib, 2, &value, &len, NULL, 0);
    value = value / 1024;
#endif

    total_mem = value * 1.e-6;
    return retcode;
}

int get_node_free_mem(double &free_mem)
{
    int retcode = 1;
    // value in KB unit
    long value = 0L;

#if defined(__linux__)
    char line[1024];
    FILE *fp = NULL;

    if ((fp = fopen("/proc/meminfo", "r")) != NULL)
    {
        while (fgets(line, sizeof(line), fp))
        {
            if (sscanf(line, "MemAvailable: %d kB", &value) == 1)
            {
                retcode = 0;
                break;
            }
        }
        fclose(fp);
    }
#endif

    free_mem = value * 1.e-6;
    return retcode;
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */

