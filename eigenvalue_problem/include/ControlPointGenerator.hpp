#ifndef CONTROLPOINTGENERATOR_HPP
#define CONTROLPOINTGENERATOR_HPP

include "BSplineBasis.hpp"
include "Eigne/Dense"

class ControlPointGenerator
{
    public:
        ControlPointGenerator() {}
        ~ControlPointGenerator() {}

        std::vector<double> GenerateControlPoints1D(
            const BSplineBasis * const &basis,
            const double &Pmax, const double &Pmin);
        
        std::vector<double> GenerateControlPoints2D(
            const BSplineBasis * const &basis1,
            const BSplineBasis * const &basis2,
            const double &P1max, const double &P1min,
            const double &P2max, const double &P2min);
}

#endif