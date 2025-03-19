#include "RefElement.hpp"

std::vector<double> Element::GenerateBasis1D(const BernsteinBasis * const &bern,
    const std::vector<double> &extraction,
    const QuadraturePoint * const &quad)
{
    std::vector<double> basis{};
    const int p = bern->GetDegree();

    const std::vector<double> BB = bern->GetBernsteinBasis(quad);

    for (int ii = 0; ii < nqp; ++ii)
    {
        std::vector<double> temp(p+1, 0.0);
        for (int jj = 0; jj < p+1; ++jj)
            for (int kk = 0; kk < p+1; ++kk)
                temp[jj] += extraction[jj*(p+1)+kk] * BB[ii*(p+1)+kk];
        basis.push_back(temp); 
    }

    return basis;
}

std::vector<double> Element::GenerateBasisDerivative1D(const BernsteinBasis * const &bern,
    const std::vector<double> &extraction,
    const QuadraturePoint * const &quad)
{
    std::vector<double> dbasis{};
    const int p = bern->GetDegree();

    const std::vector<double> dBB = bern->GetBernsteinBasisDerivative(quad);
    
    for (int ii = 0; ii < nqp; ++ii)
    {
        std::vector<double> temp(p+1, 0.0);
        for (int jj = 0; jj < p+1; ++jj)
            for (int kk = 0; kk < p+1; ++kk)
                temp[jj] += extraction[jj*(p+1)+kk] * dBB[ii*(p+1)+kk];
        dbasis.push_back(temp); 
    }

    return dbasis;
}