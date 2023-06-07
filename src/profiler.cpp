#include "profiler.h"

#include <ctime>
#include <string>
#include <cstring>
#include <iostream>
#include <omp.h>

double cpu_time_from_clocks_diff(const std::clock_t& ct_start,
                                 const std::clock_t& ct_end)
{
    return double(ct_end - ct_start) / CLOCKS_PER_SEC;
}

std::map<std::string, Profiler::Timer> Profiler::sd_map_timer;
std::map<std::string, int> Profiler::sd_map_level;
std::map<std::string, std::string> Profiler::sd_map_note;
std::vector<std::string> Profiler::sd_order;

void Profiler::Timer::start() noexcept
{
    if(is_on())
        stop();
    ncalls++;
    clock_start = clock();
    wt_start = omp_get_wtime();
    cpu_time_last = 0.0;
    wall_time_last = 0.0;
    // printf("start: %zu %zu %f\n", ncalls, clock_start, wt_start);
}

void Profiler::Timer::stop() noexcept
{
    if(!is_on()) return;
    cpu_time_accu += (cpu_time_last = cpu_time_from_clocks_diff(clock_start, clock()));
    wall_time_accu += (wall_time_last = omp_get_wtime() - wt_start);
    // printf("stop: %f %f %f\n", wt_start, wall_time, cpu_time);
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
    sd_map_timer.at(tname).start();
}

void Profiler::stop(const char *tname) noexcept
{
    if(omp_get_thread_num()!=0) return;
    if (sd_map_timer.count(tname))
    {
        sd_map_timer.at(tname).stop();
    }
    else
    {
        printf("Warning!!! Timer %s not found, profiling is very likely wrong!\n", tname);
    }
}

double Profiler::get_cpu_time_last(const char *tname) noexcept
{
    std::string sname = string(tname);
    if (sd_map_timer.count(sname))
        return sd_map_timer.at(sname).get_cpu_time_last();
    return -1.0;
}

double Profiler::get_wall_time_last(const char *tname) noexcept
{
    std::string sname = string(tname);
    if (sd_map_timer.count(sname))
        return sd_map_timer.at(sname).get_wall_time_last();
    return -1.0;
}

static std::string banner(char c, int n)
{
    std::string s = "";
    while(n--) s += c;
    return s;
}

void Profiler::display(int verbose) noexcept
{
    printf("%-49s %-12s %-18s %-18s\n", "Entry", "#calls", "CPU time (s)", "Wall time (s)");
    printf("%100s\n", banner('-', 100).c_str());
    for (auto &tname: sd_order)
    {
        // decide level indentation
        std::string s = "";
        int i = sd_map_level[tname];
        while(i--) s += "  ";
        if (verbose > 0 && i > verbose)
            continue;

        const auto name = s + sd_map_note[tname];
        const auto& t = sd_map_timer[tname];
        char cstr_cputime[100];
        char cstr_walltime[100];
        sprintf(cstr_cputime, "%.4f", t.get_cpu_time());
        sprintf(cstr_walltime, "%.4f", t.get_wall_time());
        // skip timers that have not been called
        if(!t.get_ncalls()) continue;
        printf("%-49s %-12zu %-18s %-18s\n", name.c_str(), t.get_ncalls(),
               (s + cstr_cputime).c_str(), (s + cstr_walltime).c_str());
    }
}
