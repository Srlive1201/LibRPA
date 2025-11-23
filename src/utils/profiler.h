/*!
 @file profiler.h
 @brief Utilities to profile the program
 */
#ifndef PROFILER_H
#define PROFILER_H
#include <vector>
#include <ctime>
#include <string>
#include <map>

double cpu_time_from_clocks_diff(const std::clock_t& ct_start,
                                 const std::clock_t& ct_end);

std::string get_timestamp();

//! A simple profiler object to record timing of code snippet runs in the program.
class Profiler
{
private:
    //! Class to track timing of a particular part of code
    class Timer
    {
    private:
        //! the number of timer calls
        size_t ncalls;
        //! clock when the timer is started
        std::clock_t clock_start;
        //! wall time when the timer is started
        double wt_start;
        //! accumulated cpu time
        double cpu_time_accu;
        //! accumulated wall time, i.e. elapsed time
        double wall_time_accu;
        //! cpu time during last call
        double cpu_time_last;
        //! wall time during last call
        double wall_time_last;

        // private functions
    public:
        Timer(): ncalls(0), clock_start(0), wt_start(0), cpu_time_accu(0), wall_time_accu(0), cpu_time_last(0), wall_time_last(0) {}
        //! start the timer
        void start() noexcept;
        //! stop the timer and record the timing
        void stop() noexcept;
        bool is_on() const { return clock_start != 0; };
        size_t get_ncalls() const { return ncalls; };
        double get_cpu_time() const { return cpu_time_accu; };
        double get_wall_time() const { return wall_time_accu; };
        double get_cpu_time_last() const { return cpu_time_last; };
        double get_wall_time_last() const { return wall_time_last; };
    };
    //! Container of Timer objects
    static std::map<std::string, Timer> sd_map_timer;
    //! Level of each timer to account for hierarchy
    static std::map<std::string, int> sd_map_level;
    //! Explanatory note of the timer
    static std::map<std::string, std::string> sd_map_note;
    //! Order of timers
    static std::vector<std::string> sd_order;

public:
    Profiler() = delete;
    ~Profiler() = delete;
    Profiler(const Profiler&) = delete;
    Profiler(Profiler&&) = delete;

    //! Add a timer
    static void add(const char *tname, const char *tnote = "", int level = -1) noexcept;
    //! Start a timer. If the timer is not added before, add it.
    static void start(const char *tname, const char *tnote = "", int level = -1) noexcept;
    //! Stop a timer and record the timing
    static void stop(const char *tname) noexcept;
    //! Alias to stop
    static void cease(const char *tname) noexcept;
    //! Get cpu time of last call of timer
    static double get_cpu_time_last(const char *tname) noexcept;
    //! Get wall time of last call of timer
    static double get_wall_time_last(const char *tname) noexcept;
    //! Display the current profiling result
    static void display(int verbose = 0) noexcept;
    //! Get the number of created timers
    static int get_num_timers() noexcept { return sd_order.size(); };
};

#endif
