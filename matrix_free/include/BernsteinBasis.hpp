#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

#include <vector>
#include <iostream>
#include "QuadraturePoint.hpp"

class BernsteinBasis {
    public:
        BernsteinBasis(const int &p) : p(p) {}
        ~BernsteinBasis(){}

        void GenerateBernsteinBasis(const QuadraturePoint * const &qp)
        {
            B.clear();
            dB.clear();
            const int nqp = qp->GetNumQuadraturePoint();
            const std::vector<double> &xi = qp->GetQuadraturePoint();
            for (int i = 0; i < nqp; ++i) {
                std::vector<double> tempB = GetBernsteinBasisSingleQP(xi[i]);
                std::vector<double> tempdB = GetBernsteinBasisDerivativeSingleQP(xi[i]);
                B.insert(B.end(), tempB.begin(), tempB.end());
                dB.insert(dB.end(), tempdB.begin(), tempdB.end());
            }
        }

        std::vector<double> GetBernsteinBasisSingleQP(const double &xi) const;

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const double &xi) const;

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const double &xi, const double &h) const;

        int GetDegree() const { return p; }

        const std::vector<double>& GetBernsteinBasis() const { return B; }
        const std::vector<double>& GetBernsteinBasisDerivative() const { return dB; }

    private:
        const int p;
        std::vector<double> B{};
        std::vector<double> dB{};
};

#endif