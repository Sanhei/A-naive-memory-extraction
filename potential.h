#ifndef HISTOGRAM
#define HISTOGRAM
#include <vector>
#include <math.h>

class histogram_constructor
{
//For delete class, release the memory;
//class A ...
//A* a;
//delete [] a
private:
public:
        std::vector<double> hist_plot;
        std::vector<double> hist_accum;
        //For error analyse
        std::vector<double>  free_energy_shift;

        std::vector<double> free_energy;
        std::vector<double> gradient_plot;
        std::vector<double> gradient_U;

        double temperature;
        double hist_minimum_edge;
        double hist_maximum_edge;
        int histogram_discretize;
        int hist_nbins;

        int midpoint_index;
       /*
        * The algorithm is here.
        * 1. Get the maximum and minimum of the trajectory.
        * 2. Distcretize 1 into "distcretize", e.g. if "distcretize = 100";
        *    then the bins size is 1/100.
        * 3. Accumlate it into the bins.
        */
        histogram_constructor();
        void SetTemperature(double mtemperature);
        void histogram_build(int discretize, std::vector<double>& trajectory);
        void gradient();//This function didn't join in but in case;
        void derivative(std::vector<double>& x, std::vector<double>& y, std::vector<double>& gradient);
        //The following function for analyse the sysmetic error;
        void boundary_influence(double down, double up);
        void shift_free_energy(double shift);

};

#endif
