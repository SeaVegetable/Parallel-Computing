#ifndef ELEMENTFEM_HPP
#define ELEMENTFEM_HPP

#include <vector>
#include <iostream>

class ElementFEM
{
    public:
        ElementFEM(const double &p, const double &q) : p(p), q(q) {}
        ~ElementFEM(){}

        int GetNumLocalBasis() const { return (p+1)*(q+1); }

        int GetNumLocalBasis1D(const int &dim) const
        {
            return (dim == 0) ? p+1 : q+1;
        }

        void GenerateRefElementSingleQP(const double &xi, const double &eta,
            std::vector<double> &N,
            std::vector<double> &dN_dx, std::vector<double> &dN_dy) const;
        
    private:
        const double p;
        const double q;

        std::vector<double> GenerateBasis1DSingleQP(const double &xi, const int &pp) const;

        std::vector<double> GenerateBasisDerivative1DSingleQP(const double &xi, const int &pp) const;
};

#endif