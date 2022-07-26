#include "matplotlibcpp.h"
#include <vector>
#include <cmath>
namespace plt = matplotlibcpp;
int main(){
        std::vector<double> y = {1, 2, 3, 4};
        std::vector<double> x = {3, 4, 5, 6};
        plt::figure();
        plt::plot(x, y, "bo-", {{"label", "data trend"}});
        plt::show();
        return 0;
}
