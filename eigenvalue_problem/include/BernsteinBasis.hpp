#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

#include <vector>

class BernsteinBasis {
    public:
        std::vector<double> GetBernsteinBasis1D(const int &p, const double &xi);

        std::vector<double> GetBernsteinBasis1DDerivative(const int &p, const double &xi);
};

#endif