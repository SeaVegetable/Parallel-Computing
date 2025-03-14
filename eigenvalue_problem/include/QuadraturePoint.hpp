#ifndef QUADRATUREPOINT_HPP
#define QUADRATUREPOINT_HPP

#include <vector>

class QuadraturePoint
{
    public:
        QuadraturePoint(const int &nqp, const int &left, const int &right);
        ~QuadraturePoint() {}

        std::vector<double> GetQuadraturePoint() const { return qp; }
        std::vector<double> GetWeight() const { return w; }

    private:
        std::vector<double> qp;
        std::vector<double> w;
};

#endif