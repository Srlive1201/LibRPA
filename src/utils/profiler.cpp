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

void Profiler::add(const std::string &tname, const std::string &tnote) noexcept
{
    auto new_timer = std::make_shared<Timer>(tname, tnote);

    if (!root)
    {
        root = new_timer;
    }
    else
    {
        if (current)
        {
            if (current->child)
            {
                auto sibling = current->child;
                while (sibling->next) sibling = sibling->next;
                sibling->next = new_timer;
                new_timer->prev = sibling;
            }
            else
            {
                current->child = new_timer;
            }
            new_timer->parent = current;
        }
    }

    current = new_timer;
}

void Profiler::start(const std::string &tname, const std::string &tnote) noexcept
{
    if (omp_get_thread_num() != 0) return;
    auto timer = find_timer_in_hierarchy(tname);
    // create a new timer if it is not found, otherwise we move to that timer and start
    if (!timer)
    {
        add(tname, tnote);
    }
    else
    {
        current = timer;
    }
#ifdef LIBRPA_VERBOSE
    double free_mem_gb;
    get_node_free_mem(free_mem_gb);
    global::ofs_myid << get_timestamp() <<" Timer start: " << tname << ". "
                           << "Free memory on node [GB]: " << free_mem_gb << "\n";
    std::flush(librpa_int::global::ofs_myid);
#endif
    current->start();
}

void Profiler::stop(const std::string &tname) noexcept
{
    if (omp_get_thread_num() != 0) return;
    if (current)
    {
        // Check if the current timer matches the given timer name
        if (current->name == tname)
        {
#ifdef LIBRPA_VERBOSE
            double free_mem_gb;
            get_node_free_mem(free_mem_gb);
            global::ofs_myid << get_timestamp() << " Timer stop:  " << tname << ". "
                             << "Free memory on node [GB]: " << free_mem_gb << "\n";
#endif
            current->stop();
            current = current->parent;
        }
        else
        {
             global::ofs_myid << "Warning: Attempting to stop timer '" << tname
                              << "' but current active timer is '" << current->name << "'" << std::endl;
        }
    }
    else
    {
        global::ofs_myid << "Warning: No timer is currently active" << std::endl;
    }
}

std::shared_ptr<Profiler::Timer> Profiler::find_timer_in_hierarchy(const std::string &tname)
{
    return search_timer_in_hierarchy(current, tname);
}

std::shared_ptr<Profiler::Timer> Profiler::search_timer_in_hierarchy(std::shared_ptr<Timer> timer,
                                                                     const std::string &tname)
{
    if (!timer) return nullptr;
    if (timer->name == tname) return timer;

    // Recursively check child timers
    if (timer->child)
    {
        auto found_timer = search_timer_in_hierarchy(timer->child, tname);
        if (found_timer)
        {
            return found_timer;
        }
    }

    // Check the next sibling timer
    if (timer->next)
    {
        return search_timer_in_hierarchy(timer->next, tname);
    }

    return nullptr;
}

double Profiler::get_cpu_time_last(const std::string &tname) noexcept
{
    auto timer = this->find_timer_in_hierarchy(tname);
    if (timer)
        return timer->get_cpu_time_last();
    return -1.0;
}

double Profiler::get_wall_time_last(const std::string &tname) noexcept
{
    auto timer = this->find_timer_in_hierarchy(tname);
    if (timer)
        return timer->get_wall_time_last();
    return 0.0;
}

std::string Profiler::get_profile_string_of_timer(std::shared_ptr<Profiler::Timer> timer, int level)
{
    std::ostringstream ss;
    // std::string indent(2 * level, ' ');
    std::string indent(level, ' ');
    const auto note = indent + (timer->note == "" ? timer->name : timer->note);
    std::ostringstream cstr_cputime, cstr_walltime;
    cstr_cputime << std::fixed << std::setprecision(4) << timer->get_cpu_time();
    cstr_walltime << std::fixed << std::setprecision(4) << timer->get_wall_time();

    // Print self
    ss << std::left;
    ss << std::setw(49) << note << " " << std::setw(12) << timer->ncalls << " "
        << std::setw(18) << (indent + cstr_cputime.str()) << " "
        << std::setw(18) << (indent + cstr_walltime.str()) << "\n";
    // Print child, then sibling
    if (timer->child)
        ss << get_profile_string_of_timer(timer->child, level + 1);
    if (timer->next)
        ss << get_profile_string_of_timer(timer->next, level);
    return ss.str();
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
    output << get_profile_string_of_timer(root, 0);

    return output.str();
}

namespace global
{

Profiler profiler;

}

}
