/*!
 @file profiler.h
 @brief Utilities to profile the program
 */
#ifndef PROFILER_H
#define PROFILER_H
#include <vector>
#include <ctime>
#include <string>

double cpu_time_from_clocks_diff(const std::clock_t& ct_start,
                                 const std::clock_t& ct_end);

//! A simple profiler object to record timing of code snippet runs in the program.
class Profiler
{
    private:
        //! Private class to track timing of a particular part of code
        /*!
          @note The overhead of each start/stop call is about 1e-7 s.
         */
        class Timer
        {
            private:
                // private fields
                //! the name of timer as an identifier
                std::string name;
                //! the number of timer calls
                size_t ncalls;
                //! clock when the timer is started
                std::clock_t clock_start;
                //! wall time when the timer is started
                double wt_start;
                //! accumulated cpu time
                double cpu_time;
                //! accumulated wall time, i.e. elapsed time
                double wall_time;
                // private functions
                //! check if the timer is started
                bool is_on() { return clock_start != 0; };
            public:
                Timer(const char *tn) { name = tn;
                                        ncalls = 0;
                                        clock_start = 0;
                                        wt_start = 0.;
                                        cpu_time = wall_time = 0.; };
                //! start the timer
                void start();
                //! stop the timer and record the timing
                void stop();
                std::string get_name() { return name; };
                size_t get_ncalls() { return ncalls; };
                double get_cpu_time() { return cpu_time; };
                double get_wall_time() { return wall_time; };
        };
        //! Container of Timer objects
        std::vector<Timer> timers;
        //! Level of each timer to account for hierachy
        std::vector<int> timer_levels;
        //! Explanatory note of the timer
        std::vector<std::string> timer_notes;
        //! Search the timer with requested timer name
        Timer * find_timer(const char *tname);
    public:
        Profiler(): timers(), timer_levels(), timer_notes() {};
        //! Add a timer
        void add(int level, const char *tname, const char *tnote = "");
        //! Start a timer
        void start(const char *tname);
        //! Stop a timer and record the timing
        void stop(const char *tname);
        //! Display the current profiling result
        void display();
        //! Get the number of created timers
        int get_num_timers() { return timers.size(); };
};

//! A global profiler object for convenience
extern Profiler prof;
#endif
