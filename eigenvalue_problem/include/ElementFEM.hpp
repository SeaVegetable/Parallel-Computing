#ifndef ELEMENTFEM_HPP
#define ELEMENTFEM_HPP

class ElementFEM
{
    public:
        ElementFEM(const double &p, const double &q) : p(p), q(q) {}
        ~ElementFEM(){}

        void SetElement(const double &hx, const double &hy)
        {
            this->hx = hx;
            this->hy = hy;
        }

        int GetNumLocalBasis() const { return (p+1)*(q+1); }

        int GetNumLocalBasis1D(const int &dim) const
        {
            return (dim == 0) ? p+1 : q+1;
        }

        void GenerateElementSingleQP(const double &xi, const double &eta,
            const std::vector<double> &eCP, 
            std::vector<double> &R, std::vector<double> &dR_dx, std::vector<double> &dR_dy,
            double &x, double &y,
            double &jacobian) const;
        
    private:
        const double p;
        const double q;

        double hx;
        double hy;

        std::vector<double> GenerateBasis1DSingleQP(const double &xi, const int &pp);

        std::vector<double> GenerateBasisDerivative1DSingleQP(const double &xi, const int &pp);
}

#endif