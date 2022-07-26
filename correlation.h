#ifndef CORRELATION
#define CORRELATION
#include <vector>
#include <math.h>
#include <fftw3.h>

class correlation_calculation
{
public:
     int N;// Total trajectory size-2;
     fftw_complex *velocity_;
     fftw_complex *Cor;
     std::vector<double> Correlation;
     correlation_calculation();
     void SetTime_interval(double time_interval_);
     void velocity(std::vector<double> traj);
     void fft(fftw_complex *in, fftw_complex *out);
     void ifft(fftw_complex *in, fftw_complex *out);
     void Convolution(fftw_complex *x, fftw_complex *y);
     //void Convolutiongux(fftw_complex *x, fftw_complex *y);

 private:
     double time_interval;
};
#endif
