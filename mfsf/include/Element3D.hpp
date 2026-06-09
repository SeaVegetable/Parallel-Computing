#ifndef ELEMENT3D_HPP
#define ELEMENT3D_HPP

#include <algorithm>
#include "RefElement.hpp"

class Element3D
{
    public:
        Element3D(const int &p, const int &q, const int &r) : p(p), q(q), r(r) {}
        ~Element3D(){}

        void SetElement(const std::vector<double> &extraction1,
            const std::vector<double> &extraction2,
            const std::vector<double> &extraction3,
            const double &hx, const double &hy, const double &hz)
        {
            this->extraction1 = extraction1;
            this->extraction2 = extraction2;
            this->extraction3 = extraction3;
            this->hx = hx;
            this->hy = hy;
            this->hz = hz;
        }

        int GetNumLocalBasis() const { return (p+1)*(q+1)*(r+1); }

        int GetNumLocalBasis1D(const int &dim) const
        {
            if(dim == 0)
                return p+1;
            else if (dim == 1)
                return q+1;
            else
                return r+1;
        }

        void GenerateElementSingleQP(const double &xi, const double &eta,
            const double &zeta,
            const std::vector<double> &eCP, 
            std::vector<double> &N,
            std::vector<double> &dN_dx,
            std::vector<double> &dN_dy,
            std::vector<double> &dN_dz,
            double &x, double &y, double &z,
            double &jacobian,
            std::vector<double> &invJac) const;
    
    private:
        const int p;
        const int q;
        const int r;

        std::vector<double> extraction1{};
        std::vector<double> extraction2{};
        std::vector<double> extraction3{};

        double hx;
        double hy;
        double hz;
};

#endif
