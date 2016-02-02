
#include "spline.h"
#include <cmath>
#include <vector>

/// Make a linearly spaced vector of points
template <typename T> std::vector<T> linspace(T xmin, T xmax, std::size_t n) {
    std::vector<T> x(n, 0.0);
    for ( std::size_t i = 0;  i < n; ++i) {
        x[i] = (xmax-xmin)/(n-1)*i+xmin;
    }
    return x;
}

int main()
{
    std::vector<double> x = linspace(0.0, 10.0, 10), y(10), x2 = linspace(0.0, 10.0, 200);
    for (int i = 0; i < x.size(); ++i){ y[i] = sin(x[i]); }
    for (int i = 0; i < y.size(); ++i){ std::cout << y[i] << std::endl; }
    Spline<double, double> spl(x, y);
    std::vector<double> y2 = spl.interpolate(x2);
    for (int i = 0; i < y2.size(); ++i){ std::cout << y2[i] << std::endl; }
}