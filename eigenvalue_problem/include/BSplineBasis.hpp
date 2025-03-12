#ifndef BSPLINEBASIS_HPP
#define BSPLINEBASIS_HPP

// Ref : Les Piegl and Wayne Tiller, The NURBS Book, 2nd Edition, Springer, 1997
// FindSpan and BasisFuns are based on Algorithm A2.1 and A2.2
// DerBasisFuns is based on Algorithm A2.3

#include <vector>

class BSplineBasis
{
    public:
        BSplineBasis(const int &in_p, const std::vector<double> &in_S)
            : p(in_p), S(in_S), nFunc(static_cast<int>(S.size()) - p - 1) {}
        
        ~BSplineBasis();

        int FindSpan(double &u);

        std::vector<double> BasisFuns(double &u, int &i);

        std::vector<double> DerBasisFuns(double &u, int &i, int &n);

    private:
        const int p;
        const std::vector<double> S;
        const int nFunc;
};

#endif