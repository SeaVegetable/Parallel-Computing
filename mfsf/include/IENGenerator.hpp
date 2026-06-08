#ifndef IENGENERATOR_HPP
#define IENGENERATOR_HPP

#include "BSplineBasis.hpp"

class IENGenerator
{
    public:
        std::vector<int> GenerateIEN1D(const int &nElemX, const int &p);
        std::vector<int> GenerateIEN1D(const BSplineBasis * const &basis);
        std::vector<int> GenerateIEN2D(const int &nElemX, const int &nElemY);
        std::vector<int> GenerateIEN2D(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2);
        std::vector<int> GenerateIEN3D(const int &nElemX, const int &nElemY, const int &nElemZ);
        std::vector<int> GenerateIEN3D(const BSplineBasis * const &basis1,
                                       const BSplineBasis * const &basis2,
                                       const BSplineBasis * const &basis3);
};

#endif