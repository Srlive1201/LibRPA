#include "profiler.h"

#include <ctime>
#include <string>
#include <cstring>
// #include <iostream>
#include <omp.h>
#include <chrono>
#include <sstream>
#include <iomanip>

#include "../io/global_io.h"
#ifdef LIBRPA_VERBOSE
#include "utils_mem.h"
#endif

namespace librpa_int {

double cpu_time_from_clocks_diff(const std::clock_t& ct_start,
                                 const std::clock_t& ct_end)
{
    return double(ct_end - ct_start) / CLOCKS_PER_SEC;
}

std::string get_timestamp()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    std::stringstream ss;
    ss << "[" <<std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S")
       << '.' << std::setfill('0') << std::setw(3) << milliseconds.count() << "]";
    return ss.str();
}

void Profiler::Timer::start() noexcept
{
    if(is_on())
        stop();
    ncalls++;
    clock_start = clock();
    wt_start = omp_get_wtime();
    cpu_time_last = 0.0;
    wall_time_last = 0.0;
    // librpa_int::global::lib_printf("start: %zu %zu %f\n", ncalls, clock_start, wt_start);
}

void Profiler::Timer::stop() noexcept
{
    if(!is_on()) return;
    cpu_time_accu += (cpu_time_last = cpu_time_from_clocks_diff(clock_start, clock()));
    wall_time_accu += (wall_time_last = omp_get_wtime() - wt_start);
    // librpa_int::global::lib_printf("stop: %f %f %f\n", wt_start, wall_time, cpu_time);
    wt_start = 0.;
    clock_start = 0;
}

void Profiler::add(const char *tname, const char *tnote, int level) noexcept
{
    if (sd_map_timer.count(tname) == 0)
    {
        // when level is negative, set the level according to status of previous timers
        if (level < 0)
        {
            level = 0;
            for (auto it = sd_order.crbegin(); it != sd_order.crend(); it++)
            {
                const auto t = sd_map_timer[*it];
                if (t.is_on())
                {
                    level = sd_map_level[*it] + 1;
                    break;
                }
            }
        }
        sd_map_timer[tname] = Timer();
        sd_map_level[tname] = level;
        if (strcmp(tnote, "") == 0)
            sd_map_note[tname] = tname;
        else
            sd_map_note[tname] = tnote;
        sd_order.push_back(tname);
    }
}

void Profiler::start(const char *tname, const char *tnote, int level) noexcept
{
    if(omp_get_thread_num()!=0) return;
    if (sd_map_timer.count(tname) == 0)
    {
        add(tname, tnote, level);
    }
#ifdef LIBRPA_VERBOSE
    double free_mem_gb;
    get_node_free_mem(free_mem_gb);
    global::ofs_myid << get_timestamp() <<" Timer start: " << tname << ". "
                           << "Free memory on node [GB]: " << free_mem_gb << "\n";
    std::flush(librpa_int::global::ofs_myid);
#endif
    sd_map_timer.at(tname).start();
}

void Profiler::stop(const char *tname) noexcept
{
    if(omp_get_thread_num()!=0) return;
    if (sd_map_timer.count(tname))
    {
#ifdef LIBRPA_VERBOSE
        double free_mem_gb;
        get_node_free_mem(free_mem_gb);
        global::ofs_myid << get_timestamp() << " Timer stop:  " << tname << ". "
                               << "Free memory on node [GB]: " << free_mem_gb << "\n";
#endif
        sd_map_timer.at(tname).stop();
    }
    else
    {
        librpa_int::global::lib_printf("Warning!!! Timer %s not found, profiling is very likely wrong!\n", tname);
    }
}

double Profiler::get_cpu_time_last(const char *tname) noexcept
{
    std::string sname(tname);
    if (sd_map_timer.count(sname))
        return sd_map_timer.at(sname).get_cpu_time_last();
    return -1.0;
}

double Profiler::get_wall_time_last(const char *tname) noexcept
{
    std::string sname(tname);
    if (sd_map_timer.count(sname))
        return sd_map_timer.at(sname).get_wall_time_last();
    return 0.0;
}

static std::string banner(char c, int n)
{
    std::string s = "";
    while(n--) s += c;
    return s;
}

void Profiler::display(int verbose) noexcept
{
    global::lib_printf("%s", this->get_profile_string(verbose).c_str());
}

std::string Profiler::get_profile_string(int verbose) noexcept
{
    std::ostringstream output;
    output << std::left;

    output << std::setw(49) << "Entry" << " " << std::setw(12) << "#calls" << " "
        << std::setw(18) << "CPU time (s)" << " " << std::setw(18) << "Wall time (s)" << "\n";
    output << banner('-', 100) << "\n";

    for (auto &tname: sd_order)
    {
        int i = sd_map_level[tname];
        if (verbose > 0 && i > verbose)
            continue;

        // Decide level indentation
        std::string indent = "";
        while (i--) indent += "  ";

        const auto name = indent + sd_map_note[tname];
        const auto& t = sd_map_timer[tname];
        std::ostringstream cstr_cputime, cstr_walltime;

        cstr_cputime << std::fixed << std::setprecision(4) << t.get_cpu_time();
        cstr_walltime << std::fixed << std::setprecision(4) << t.get_wall_time();

        // Skip timers that have not been called
        if (!t.get_ncalls()) continue;

        output << std::setw(49) << name << " " << std::setw(12) << t.get_ncalls() << " "
            << std::setw(18) << (indent + cstr_cputime.str()) << " "
            << std::setw(18) << (indent + cstr_walltime.str()) << "\n";
    }

    return output.str();
}

namespace global
{

Profiler profiler;

}

}
