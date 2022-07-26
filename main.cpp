#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include "potential.h"
#include "correlation.h"
#include <cmath>
#include "spline.h"
#include <fftw3.h>
#include "figure_plot.h"
#include <typeinfo>
#include "memory.h"
//For spline
#include <boost/math/interpolators/makima.hpp>
void analytical_potential(std::vector<double> &traj, std::vector<double> &potential_a)
{
    double depth = 3.0;
    for(unsigned int i; i<traj.size(); i++){
        potential_a.push_back(depth*(4*pow(traj[i], 3) - 4*traj[i]));
    }
}
int main(int argc, char **argv)
{
    double time_interval = 0.001;
    double temperature   = 1.0;
    std::ifstream inputfile;
    inputfile.open("trajectory.txt");
    if(!inputfile)
    {
        std::cout<<"unable to openfile";
        exit(1);
    }
    double x, y;
    //Since in the trajectory one column is time, another is position;
    std::vector<double> trajectory;
    std::cout<<"Reading documents"<<std::endl;
    while(inputfile>>x>>y)
    {
        trajectory.push_back(y);
    }



    figure_plot Figure;//This class we called matplotlib.pyplot function from python
    /*
     * Histogram bin parameter, the histogram will depends on
     * your data. this discretize means each bin size.
     * For example. [1.4, 1.2, 1.3, 1.4, 1.5, 1.2],
     * the min is 1.2, and if you choose distcretize=10,
     * then there will be [1.5-1.2]*10= 3 bins,
     * if distcretize = 100, there will be 30 bins etc.
     */
    int distcretize = 50;
    histogram_constructor Hist_Obj;
    Hist_Obj.SetTemperature(temperature);
    Hist_Obj.histogram_build(distcretize, trajectory);
    Hist_Obj.boundary_influence(-1.5, 1.6);
    //Here we set the boundary is -1.5 and 1.5;
    std::vector<double> Free_energy;
    for(unsigned int i; i<Hist_Obj.hist_accum.size(); i++)
    {
        Free_energy.push_back(-temperature*log(Hist_Obj.hist_accum[i]));
    }
    std::cout<<"Histogram end"<<std::endl;
    std::cout<<"The boundary influnce"<<std::endl;
    /*
     * In this part we will integrate over the boundary.
     * This case, our well position is 1 and -1, so we set
     * the boundary is -1.5 and 1.5.
     * The method to study is to integrate the area beyond
     * this part([-1.5, 1.5]) over the whole integrate.
     * So we can see the rate of this part in the trajectory.
     */





     //For theoretical plot
    std::vector<double> theoretical;
    std::cout<<"Size is "<<Hist_Obj.hist_plot.size()<<std::endl;
    std::vector<double> plot_edge;

    double shift; //for error analyse, make the simulation result and theoretical overlap
    for(unsigned int i=0; i < Hist_Obj.hist_plot.size(); i++)
    {
        plot_edge.push_back(Hist_Obj.hist_plot[i]);
        theoretical.push_back(3.0*(pow(Hist_Obj.hist_plot[i],4)-2*pow(Hist_Obj.hist_plot[i], 2))+3);
        if(i == Hist_Obj.midpoint_index)
                shift = theoretical[i]-Hist_Obj.free_energy[i];
    }
    // Do the shift and see the result
    Hist_Obj.shift_free_energy(shift);



//Get the function of Free energy( spline )
    /* Two spline mthods:
     * 1. Cubic spline, use 3 order polynomiers ax^3+bx^2+cx+d;
     *    tk::spline, need to call spline.h file.
     * 2. Akima spline, the local suddenly drop or raise doesn't
     *    change the other distinct, which means, the function has
     *    less vibration;
     */
    tk::spline Potential;
    Potential.set_points(Hist_Obj.hist_plot, Free_energy, tk::spline::cspline);
    //using boost::math::interpolators::makima;
    //auto Potential = makima<std::vector<double>>(std::move(Hist_Obj.hist_plot), std::move(Hist_Obj.free_energy));
    //Spline function saved in Potential;
       //For spline plot
    int splinesize = int((Hist_Obj.hist_maximum_edge - Hist_Obj.hist_minimum_edge)/(double)0.0001);
    //std::cout<<"Spline fault"<<std::endl;
    std::vector<double> smallbins;
    std::vector<double> Energy;
    for(unsigned int i=0; i < splinesize; i++)
    {
        smallbins.push_back(Hist_Obj.hist_minimum_edge+0.0001*i);
        Energy.push_back(Potential(Hist_Obj.hist_minimum_edge+0.0001*i));
    }
    /* This will plot the potential, and smooth function
     * smallbins is the for spline function, show the detail of
     * spline result.
     */
    //Figure.figure_plot::potential_plot(plot_edge,smallbins, Energy, Free_energy, theoretical, Hist_Obj.free_energy_shift);

    //Clear the memory of shift free energy;
    Hist_Obj.free_energy_shift.clear();


    std::cout<<typeid(Potential).name()<<" is the type of Potential"<<std::endl;



    //Gradient of potentail
    std::cout<<"Gradient of potentail"<<std::endl;
    std::vector<double> gradient_u_plot;
    std::vector<double> analytical_gradient;
    std::vector<double> gradient_spline_plot;
    Hist_Obj.derivative(plot_edge, Free_energy, gradient_u_plot);

    analytical_potential(plot_edge, analytical_gradient);


    for(unsigned int i=0; i<plot_edge.size(); i++){
        gradient_spline_plot.push_back(plot_edge[i]);
    }
    //Spline;
    //using boost::math::interpolators::makima;
    //auto G_u = makima<std::vector<double>>(std::move(plot_edge), std::move(gradient_u_plot));
    tk::spline G_u;
    G_u.set_points(plot_edge, gradient_u_plot, tk::spline::cspline);

    std::cout<<"spline gradient right"<<std::endl;
    std::vector<double> gradient_u_spline;
    for(unsigned int i=0; i<gradient_spline_plot.size(); i++){
        if(gradient_spline_plot[i]>2)
                std::cout<<"This number is rong"<<std::endl;
        gradient_u_spline.push_back(G_u(gradient_spline_plot[i]));
    }
    std::cout<<"right"<<std::endl;
   // Figure.figure_plot::gradient_plot(gradient_spline_plot, gradient_u_spline, analytical_gradient);


    //Plot gradien
    //Release memory
    Free_energy.clear();
    theoretical.clear();
    plot_edge.clear();
    smallbins.clear();
    Energy.clear();
    gradient_u_plot.clear();
    analytical_gradient.clear();
    gradient_spline_plot.clear();
    gradient_u_spline.clear();

    Hist_Obj.gradient_plot.clear();
    Hist_Obj.gradient_U.clear();
    Hist_Obj.free_energy_shift.clear();










    //Correlation Calculation
    std::cout<<"Correlation of velocity_and velocity"<<std::endl;
    correlation_calculation Cor_vv;
    Cor_vv.N = trajectory.size() - 2;
    Cor_vv.SetTime_interval(time_interval);
    Cor_vv.correlation_calculation::velocity(trajectory);
    Cor_vv.correlation_calculation::fft(Cor_vv.velocity_, Cor_vv.velocity_);
    Cor_vv.correlation_calculation::Convolution(Cor_vv.velocity_, Cor_vv.velocity_);
    //Figure.figure_plot::Corvv_plot(time_interval, Cor_vv.Correlation);
    std::ofstream corvv_r;
    corvv_r.open("corvv.txt");
    for(size_t i; i<Cor_vv.Correlation.size();++i)
            corvv_r<<Cor_vv.Correlation[i]<<std::endl;
    corvv_r.close();

    Cor_vv.Correlation.clear();







    //End of the  Corvv and begin the Cor_gux
    std::cout<<"The correlation of gradient Potential and position"<<std::endl;

    correlation_calculation Cor_gux;
    Cor_gux.SetTime_interval(time_interval);
    Cor_gux.N = trajectory.size() - 2;
    //Transfer the vector to fftw_complex
    static fftw_complex *trajectory_f;
    trajectory_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(trajectory.size()-2));
    for(unsigned int i = 0; i<trajectory.size()-2; i++)
    {
        trajectory_f[i][0] = trajectory[i+1];
    }
    static fftw_complex *gradient_u;
    gradient_u = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(trajectory.size()-2));
    double depth = 3.0;
    for(unsigned int i=0; i<trajectory.size()-2; i++)
    {
        gradient_u[i][0] = G_u(trajectory[i+1]);
        //gradient_u[i][0] = depth*(4*pow(trajectory[i], 3) - 4*trajectory[i]);
        //gradient_u[i][0] = spline(trajectory[i]);
    }
    trajectory.clear();
    Cor_gux.correlation_calculation::fft(trajectory_f, trajectory_f);
    Cor_gux.correlation_calculation::fft(gradient_u, gradient_u);
    Cor_gux.correlation_calculation::Convolution(gradient_u, trajectory_f);
    //Figure.figure_plot::Corg_ux_plot(time_interval, Cor_gux.Correlation);
    std::ofstream corgux_r;
    std::cout<<"Recording Corux"<<std::endl;
    corgux_r.open("corux.txt");
    for(size_t i; i<Cor_gux.Correlation.size();++i)
            corgux_r<<Cor_gux.Correlation[i]<<std::endl;
    corgux_r.close();
    //Figure.figure_plot::Correlation_plot(time_interval, Cor_gux.Correlation, Cor_vv.Correlation);
    std::cout<<"Untill"<<std::endl;
    return 0;

}
