#ifndef IENGENERATOR_HPP
#define IENGENERATOR_HPP

#include "BSplineBasis.hpp"

class IENGenerator
{
    public:
        std::vector<int> GenerateIEN1D(const int &nElemX, const int &p);
        std::vector<int> GenerateIEN1D(const BSplineBasis * const &basis);
        std::vector<int> GenerateIEN2D(const int &nElemX, const int &nElemY, const int &p, const int &q);
        std::vector<int> GenerateIEN2D(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2);
};

#endif