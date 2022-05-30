#include <ctime>
#include <string>
#include <iostream>
#include <omp.h>
#include "profiler.h"

void Profiler::Timer::start()
{
    if(is_on())
        stop();
    ncalls++;
    clock_start = clock();
    wt_start = omp_get_wtime();
    // printf("start: %zu %zu %f\n", ncalls, clock_start, wt_start);
}

void Profiler::Timer::stop()
{
    if(!is_on()) return;
    cpu_time += double(clock() - clock_start) / CLOCKS_PER_SEC;
    wall_time += omp_get_wtime() - wt_start;
    // printf("stop: %f %f %f\n", wt_start, wall_time, cpu_time);
    wt_start = 0.;
    clock_start = 0;
}

Profiler::Timer * Profiler::find_timer(const char *tname)
{
    Profiler::Timer * pt = nullptr;
    for (auto i = 0; i != timers.size(); i++)
    {
        if(timers[i].get_name() == tname)
        {
            pt = &timers[i];
            // printf("%s\n", i.get_name().c_str());
            break;
        }
    }
    return pt;
}

void Profiler::add(int level, const char *tname, const char *tnote)
{
    timer_levels.push_back(level);
    timers.push_back(Timer(tname));
    // printf("add: %zu\n", t.get_ncalls());
    std::string note = tnote;
    if(note == "")
        timer_notes.push_back(std::string(tname));
    else
        timer_notes.push_back(std::string(tnote));
}

void Profiler::start(const char *tname)
{
    if(omp_get_thread_num()!=0) return;
    Timer * pt = find_timer(tname);
    if(pt!=nullptr)
        pt->start();
}

void Profiler::stop(const char *tname)
{
    if(omp_get_thread_num()!=0) return;
    Timer * pt = find_timer(tname);
    if(pt!=nullptr)
        pt->stop();
}

std::string banner(char c, int n)
{
    std::string s = "";
    while(n--) s += c;
    return s;
}

void Profiler::display()
{
    printf("%-45s %14s %19s %19s\n", "Entry", "#calls", "CPU time (s)", "Wall time (s)");
    printf("%100s\n", banner('-', 100).c_str());
    Timer * pt;
    for (auto t = 0; t < timers.size(); t++)
    {
        std::string s = "";
        int i = timer_levels[t];
        while(i--) s += "  ";
        s += timer_notes[t];
        pt = &timers[t];
        // skip timers that have not been called
        if(!pt->get_ncalls()) continue;
        printf("%-45s %14zu %19.4f %19.4f\n", s.c_str(), pt->get_ncalls(), pt->get_cpu_time(), pt->get_wall_time());
    }
}

Profiler prof;
