#ifndef REFELEMENT_HPP
#define REFELEMENT_HPP

#include "BernsteinBasis.hpp"

class RefElement
{
    public:
        std::vector<double> GenerateBasis(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const QuadraturePoint * const &quad);
        
        std::vector<double> GenerateBasisDerivative(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const QuadraturePoint * const &quad);
};

#endif