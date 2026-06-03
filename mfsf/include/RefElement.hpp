#ifndef MFSF_REF_ELEMENT_HPP
#define MFSF_REF_ELEMENT_HPP

#include <vector>

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
            const double &xi,
            const double &h);

        std::vector<double> GenerateBasisDerivative1DSingleQP(const std::vector<double> &dB,
            const std::vector<double> &extraction);
};

#endif
