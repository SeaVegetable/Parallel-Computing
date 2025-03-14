#ifndef NURBSEXTRACTIONGENERATOR_HPP
#define NURBSEXTRACTIONGENERATOR_HPP

// Ref : Borden et al. 2011 IJNME, 87:15-47
// GenerateExtraction1D : Algorithm 1

#include "BSplineBasis.hpp"

class NURBSExtractionGenerator
{
    public:
        std::vector<int> GenerateExtraction1D(const BSplineBasis * const &basis);

    private:
        inline int idxC(int &i, int &j, int &k, int &a)
        {
            return i + j * a + k * a * a;
        }
};

#endif