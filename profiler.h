#ifndef PROFILER_H
#define PROFILER_H
#include <vector>
#include <ctime>
#include <string>

//! A simple profiler object to record timing of code snippet runs in the program.
class Profiler
{
    private:
        //! Private class to track timing of a particular part of code
        /*!
          \note The overhead of each start/stop call is about 1e-7 s.
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
        std::vector<Timer> timers;
        //! level of each timer to account for hierachy
        std::vector<int> timer_levels;
        //! explanatory note of the timer
        std::vector<std::string> timer_notes;
        //! search the timer with requested timer name
        Timer * find_timer(const char *tname);
    public:
        Profiler(): timers(), timer_notes(), timer_levels() {};
        //! add a timer
        void add(int level, const char *tname, const char *tnote = "");
        //! start a timer
        void start(const char *tname);
        //! stop a timer and record the timing
        void stop(const char *tname);
        //! display the profiling result
        void display();
        int get_num_timers() { return timers.size(); };
};

//! a global profiler object
extern Profiler prof;
#endif
