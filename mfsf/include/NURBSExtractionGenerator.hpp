#ifndef MFSF_NURBS_EXTRACTION_GENERATOR_HPP
#define MFSF_NURBS_EXTRACTION_GENERATOR_HPP

#include <vector>

#include "BSplineBasis.hpp"

class NURBSExtractionGenerator
{
    public:
        std::vector<double> GenerateExtraction1D(const BSplineBasis * const &basis);

    private:
        int idxC(const int &i, const int &j, const int &k, const int &a) const
        {
            return (j - 1) + (i - 1) * a + (k - 1) * a * a;
        }
};

#endif
