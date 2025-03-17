#ifndef IENGENERATOR_HPP
#define IENGENERATOR_HPP

#include "BSplineBasis.hpp"

class IENGenerator
{
    public:
        std::vector<int> GenerateIEN1D(const BSplineBasis * const &basis);
        std::vector<int> GenerateIEN2D(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2);
};

#endif