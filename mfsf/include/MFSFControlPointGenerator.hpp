#ifndef MFSF_CONTROL_POINT_GENERATOR_HPP
#define MFSF_CONTROL_POINT_GENERATOR_HPP

#include <array>
#include <stdexcept>
#include <vector>

#include "BSplineBasis.hpp"

class MFSFControlPointGenerator
{
    public:
        std::vector<double> GenerateControlPoints1D(const BSplineBasis * const &basis,
            const double &min_value,
            const double &max_value) const;

        template<int Dim>
        std::vector<double> GenerateControlPoints(
            const std::array<const BSplineBasis *, Dim> &bases,
            const std::array<double, Dim> &min_corner,
            const std::array<double, Dim> &max_corner) const
        {
            std::array<std::vector<double>, Dim> cp_axes;
            std::array<int, Dim> axis_sizes{};
            int total_points = 1;

            for (int axis = 0; axis < Dim; ++axis)
            {
                if (bases[axis] == nullptr)
                {
                    throw std::invalid_argument("BSpline basis pointer cannot be null");
                }
                cp_axes[axis] = GenerateControlPoints1D(bases[axis], min_corner[axis], max_corner[axis]);
                axis_sizes[axis] = static_cast<int>(cp_axes[axis].size());
                total_points *= axis_sizes[axis];
            }

            std::vector<double> control_points(total_points * Dim, 0.0);
            std::array<int, Dim> multi_index{};
            for (int point_id = 0; point_id < total_points; ++point_id)
            {
                DecodeIndex<Dim>(point_id, axis_sizes, multi_index);
                for (int axis = 0; axis < Dim; ++axis)
                {
                    control_points[point_id * Dim + axis] = cp_axes[axis][multi_index[axis]];
                }
            }

            return control_points;
        }

    private:
        std::vector<double> SolveLinearSystem(std::vector<double> A,
            std::vector<double> b,
            const int n) const;

        template<int Dim>
        void DecodeIndex(const int flat,
            const std::array<int, Dim> &widths,
            std::array<int, Dim> &multi_index) const
        {
            int remainder = flat;
            for (int axis = 0; axis < Dim; ++axis)
            {
                multi_index[axis] = remainder % widths[axis];
                remainder /= widths[axis];
            }
        }
};

#endif
