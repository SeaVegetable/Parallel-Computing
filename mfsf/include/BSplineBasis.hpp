#ifndef MFSF_BSPLINE_BASIS_HPP
#define MFSF_BSPLINE_BASIS_HPP

#include <vector>

class BSplineBasis
{
    public:
        BSplineBasis(const int &in_p, const std::vector<double> &in_S)
            : p(in_p), S(in_S) {}

        int GetDegree() const { return p; }

        std::vector<double> GetKnotVector() const { return S; }

        int GetNumFunctions() const
        {
            return static_cast<int>(S.size()) - p - 1;
        }

        int GetNumElements() const
        {
            return GetNumFunctions() - p;
        }

        int FindSpan(const double &u) const;

        std::vector<double> BasisFuns(const double &u, const int &i) const;

        std::vector<double> DerBasisFuns(const double &u, const int &i, const int &n) const;

    private:
        const int p;
        const std::vector<double> S;
};

#endif
