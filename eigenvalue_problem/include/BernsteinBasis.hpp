#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

#include QuadraturePoint.hpp

class BernsteinBasis {
    public:
        BernsteinBasis(const int &p, const QuadraturePoint * const &quad);
        ~BernsteinBasis(){}

        std::vector<double> GetBernsteinBasisSingleQP(const int &p, const double &xi);

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const int &p, const double &xi);

        std::vector<double> GetBernsteinBasis() const { return BB; }

        std::vector<double> GetBernsteinBasisDerivative() const { return dBB; }

    private:
        std::vector<double> BB{};
        std::vector<double> dBB{};
};

#endif