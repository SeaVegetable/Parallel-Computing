#ifndef MFSF_QUADRATURE_POINT_HPP
#define MFSF_QUADRATURE_POINT_HPP

#include <vector>

class QuadraturePoint
{
    public:
        QuadraturePoint(const int &nqp, const double &left, const double &right);

        int GetNumQuadraturePoint() const { return nqp; }
        std::vector<double> GetQuadraturePoint() const { return qp; }
        std::vector<double> GetWeight() const { return w; }

    private:
        const int nqp;
        std::vector<double> qp;
        std::vector<double> w;
};

#endif
