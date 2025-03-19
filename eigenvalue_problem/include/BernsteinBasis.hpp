#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

#include QuadraturePoint.hpp

class BernsteinBasis {
    public:
        BernsteinBasis(const int &p) : p(p) {}
        ~BernsteinBasis(){}

        std::vector<double> GetBernsteinBasisSingleQP(const double &xi);

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const double &xi);

        int GetDegree() const { return p; }

        std::vector<double> GetBernsteinBasis(const QuadraturePoint * const &quad);

        std::vector<double> GetBernsteinBasisDerivative(const QuadraturePoint * const &quad);

    private:
        const int p;
};

#endif