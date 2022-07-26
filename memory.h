#ifndef MEMORY
#define MEMORY

#include <vector>
#include <math.h>

class memory
{
public:
    memory();
    std::vector<double> G;
    void memory_calculation(double time_interval, std::vector<double> corvv, std::vector<double> corux, int memory_range);
};
#endif
