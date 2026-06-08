#ifndef IDGENERATOR_HPP
#define IDGENERATOR_HPP

#include "BSplineBasis.hpp"

class IDGenerator
{
    public:
        std::vector<int> GenerateID1D(const int &nFuncX);
        std::vector<int> GenerateID1D(const BSplineBasis * const &basis);
        std::vector<int> GenerateID2D(const int &nFuncX, const int &nFuncY);
        std::vector<int> GenerateID2D(const BSplineBasis * const &basis1, const BSplineBasis * const &basis2);
        std::vector<int> GenerateID3D(const int &nFuncX, const int &nFuncY, const int &nFuncZ);
        std::vector<int> GenerateID3D(const BSplineBasis * const &basis1,
                                      const BSplineBasis * const &basis2,
                                      const BSplineBasis * const &basis3);
};

#endif