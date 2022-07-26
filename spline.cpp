#include <boost/math/interpolators/makima.hpp>
#include <vector>

int main(){
    std::vector<double> x{1, 2, 3 , 4};
    std::vector<double> y{1,4, 9, 16};
    using boost::math::interpolators::makima;
    auto spline = makima<std::vector<double>>(std::move(x), std::move(y));
// evaluate at a point:
    double z = spline(3.4);
// evaluate derivative at a point:
    double zprime = spline.prime(3.4);
    std::cout<<"spline result"<<z<<std::endl;

    std::cout<<"spline derivative result"<<zprime<<std::endl;
    return 0;
}

