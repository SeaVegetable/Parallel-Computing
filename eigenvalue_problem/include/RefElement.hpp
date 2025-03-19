#ifndef REFELEMENT_HPP
#define REFELEMENT_HPP

#include "BernsteinBasis.hpp"

class RefElement
{
    public:
        std::vector<double> GenerateBasis1D(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const QuadraturePoint * const &quad);
        
        std::vector<double> GenerateBasisDerivative1D(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const QuadraturePoint * const &quad);
};

#endif