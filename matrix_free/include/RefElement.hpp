#ifndef REFELEMENT_HPP
#define REFELEMENT_HPP

#include "BernsteinBasis.hpp"

class RefElement
{
    public:
        std::vector<double> GenerateBasis1DSingleQP(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const double &xi);
        
        std::vector<double> GenerateBasis1DSingleQP(const std::vector<double> &B,
            const std::vector<double> &extraction);
        
        std::vector<double> GenerateBasisDerivative1DSingleQP(const BernsteinBasis * const &bern,
            const std::vector<double> &extraction,
            const double &xi, const double &h);
        
        std::vector<double> GenerateBasisDerivative1DSingleQP(const std::vector<double> &dB,
            const std::vector<double> &extraction);
};

#endif