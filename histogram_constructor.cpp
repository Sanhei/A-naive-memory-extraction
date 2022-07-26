#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include "potential.h"
#include "spline.h"
histogram_constructor :: histogram_constructor(){}
void histogram_constructor::SetTemperature(double mtemperature)
{
    temperature = mtemperature;
}



void histogram_constructor::histogram_build(int discretize, std::vector<double>& trajectory)
{
    histogram_discretize = discretize;
    double max_boundary = *std::max_element(trajectory.begin(), trajectory.end());
    double min_boundary = *std::min_element(trajectory.begin(), trajectory.end());
    std::vector<int> histogram_init;
    for(unsigned int i = 0; i < trajectory.size();i++ )
    {
        histogram_init.push_back(static_cast<int>(std::floor(trajectory[i]*histogram_discretize+0.5)));
    }
    int hist_min = std::floor(min_boundary*histogram_discretize+0.5);
    int hist_max = std::floor(max_boundary*histogram_discretize+0.5);
    hist_nbins = hist_max - hist_min;
    for(unsigned int i=0; i<hist_nbins; i++)
    {
        hist_accum.push_back(0);
    }
    std::cout<<"The histogram maximum:"<<hist_max/(double)histogram_discretize<<std::endl;
    hist_maximum_edge = hist_max/(double)histogram_discretize;
    std::cout<<"The number of bins is "<<hist_nbins<<std::endl;
    std::cout<<"The bins' size is "<<1/(double)discretize<<std::endl;
    //add a pseudo point in the end, max+1, the reason because we use floor function,

    for(unsigned int i = 0; i<hist_nbins; i++)
        hist_plot.push_back(0);
    for(int i = 0; i<hist_nbins; i++)
    {
        hist_plot[i]=(hist_min+i)/(double)histogram_discretize;
        if(i==0)
        {
            std::cout<<"Test: the minmum edge of histogram is "<<hist_min/(double)histogram_discretize<<std::endl;
            hist_minimum_edge = hist_min/(double)histogram_discretize;
        }
        if(hist_min+i<0) midpoint_index = i;
        //This variable for error analysis, in the end we need to shift a bit
        //to compare the shape with theoretical result; subtract the difference
        //between the zero point.
    }
    //This gives each point probability weight.
    static double rescale =  1*histogram_discretize/(double)trajectory.size();
    int hist_index = 0;
    for(unsigned int i= 0; i < histogram_init.size(); i++)
    {
        hist_accum[histogram_init[i]-hist_min] +=  1*rescale;
    }
    histogram_init.clear();
    //consider boundary
   if(max_boundary>hist_plot.back())
    {
        hist_plot.push_back(hist_plot.back()+1/(double)histogram_discretize);
        hist_accum.push_back(hist_accum.back()/(double)2);
    }
    if(min_boundary<hist_plot[0])
    {
        hist_plot.insert(hist_plot.begin(), hist_plot[0]-1/(double)histogram_discretize );
        midpoint_index += 1;
        hist_accum.insert(hist_accum.begin(), hist_accum[0]/2 );
    }//For the max boundary


    //Test out:
    double sum = 0;
    for( unsigned int i = 0; i<hist_accum.size(); i++)
    {
        sum += hist_accum[i]/histogram_discretize;
    }



    std::cout<<"Integrate over the histogram is "<<sum<<std::endl;




    //Free energy calculation: F = -k_BTlog(p);
    for(int j = 0; j<hist_plot.size(); j++)
    {
        free_energy.push_back(0);
    }
    std::cout<<"Temperature is "<<temperature<<std::endl;
    for(int j=0; j<hist_accum.size(); j++)
    {
        if(hist_accum[j]==0)
        {
            free_energy[j] = 0.0;
        }
        else
        {
            free_energy[j] = -temperature*log(hist_accum[j]);
        }
    }
    std::cout<<"size is "<<hist_plot.size()<<std::endl;
}



void histogram_constructor::gradient()
{
//for boundary
    gradient_plot.push_back(hist_plot[0]);
    gradient_U.push_back((free_energy[1]-free_energy[0])*histogram_discretize-((free_energy[2]-free_energy[1])-(free_energy[1]-free_energy[0]))*histogram_discretize*0.5);
    std::cout<<"right"<<std::endl;
    //Use second derivative: x[0]' = x[0.5]' - x"*0.5h
    for(unsigned int i=0; i<hist_plot.size()-1; i++ )
    {
        gradient_plot.push_back((hist_plot[i+1]+hist_plot[i])/2);
        //Get mid point of each point of bins
        gradient_U.push_back((free_energy[i+1]-free_energy[i])/(hist_plot[i+1]-hist_plot[i]));
        //(x(i+1)-x(i))/stepsize
    }
    //for boundary
    gradient_plot.push_back(hist_plot.back());
    gradient_U.push_back((free_energy.end()[-2]-free_energy[1])*histogram_discretize+((free_energy.end()[-3]-free_energy.end()[-2])+(free_energy.end()[1-2]-free_energy.end()[-1]))*histogram_discretize*0.5);
}



void histogram_constructor::derivative(std::vector<double>& x, std::vector<double>& y, std::vector<double>& gradient)
{
//for boundary at beginning
    if((x[1]-x[0])==0)
        std::cout<<"Need to check, dominator equal to zero"<<std::endl;
    gradient.push_back((y[1]-y[0])/(x[1]-x[0]));
    //The idea is we use the two points beside the x[i], so which means get the mid point derivative;
    for(unsigned int i=1; i<x.size()-1; i++ )
    {
        //Get mid point of each point of bins
        gradient.push_back((y[i+1]-y[i-1])/(x[i+1]-x[i-1]));
        if((x[i+1]-x[i-1])==0)
            std::cout<<"Need to check, dominator equal to zero"<< i <<std::endl;
       //(x(i+1)-x(i))/stepsize
    }
    //For boundary of the end;
    gradient.push_back((y.end()[-1]-y.end()[-2])/(x.end()[-1]-x.end()[-2]));
    if((x.end()[-1]-x.end()[-2])==0)
        std::cout<<"Need to check, dominator equal to zero"<<"the end"<<std::endl;
}




void histogram_constructor::boundary_influence(double down, double up)
{
    double sum = 0.0;
    //Get the index for x<-1.5 and x > 1.5;
    //The sum the corresponding probability(hist_accum);
    for(unsigned int i = 0; i < hist_plot.size(); ++i)
    {
        if(hist_plot[i]<down) sum += hist_accum[i];
        if(hist_plot[i]>  up) sum += hist_accum[i];
    }
    sum /= histogram_discretize;
    std::cout<<"The probability for integrate over the boundary over the whole trajectory is "<< sum<<std::endl;
}

void histogram_constructor::shift_free_energy(double shift)
//Somehow sometime there may diverge a bit between the theoretical and simulation;
{
     for(auto& value : free_energy)
     {
        free_energy_shift.push_back(value + shift);
     }
}
