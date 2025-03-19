#ifndef BERNSTEINBASIS_HPP
#define BERNSTEINBASIS_HPP

class BernsteinBasis {
    public:
        BernsteinBasis(const int &p) : p(p) {}
        ~BernsteinBasis(){}

        std::vector<double> GetBernsteinBasisSingleQP(const double &xi);

        std::vector<double> GetBernsteinBasisDerivativeSingleQP(const double &xi, const double &h);

        int GetDegree() const { return p; }

    private:
        const int p;
};

#endif