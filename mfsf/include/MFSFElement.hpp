#ifndef MFSF_ELEMENT_HPP
#define MFSF_ELEMENT_HPP

#include <array>
#include <stdexcept>
#include <vector>

template<int Dim>
struct MFSFEvaluation
{
    int n_basis;
    std::vector<double> basis_values;
    std::vector<double> basis_gradients_param;
    std::vector<double> basis_gradients_physical;
    std::array<double, Dim> physical_point;
    std::array<double, Dim * Dim> jacobian;
    std::array<double, Dim * Dim> inv_jacobian;
    double det_jacobian;
};

template<int Dim>
struct JacobianOps;

template<>
struct JacobianOps<2>
{
    static double Determinant(const std::array<double, 4> &j)
    {
        return j[0] * j[3] - j[1] * j[2];
    }

    static std::array<double, 4> Inverse(const std::array<double, 4> &j, const double det)
    {
        std::array<double, 4> inv{};
        inv[0] = j[3] / det;
        inv[1] = -j[1] / det;
        inv[2] = -j[2] / det;
        inv[3] = j[0] / det;
        return inv;
    }
};

template<>
struct JacobianOps<3>
{
    static double Determinant(const std::array<double, 9> &j)
    {
        return j[0] * (j[4] * j[8] - j[5] * j[7])
            - j[1] * (j[3] * j[8] - j[5] * j[6])
            + j[2] * (j[3] * j[7] - j[4] * j[6]);
    }

    static std::array<double, 9> Inverse(const std::array<double, 9> &j, const double det)
    {
        std::array<double, 9> inv{};
        inv[0] = (j[4] * j[8] - j[5] * j[7]) / det;
        inv[1] = (j[2] * j[7] - j[1] * j[8]) / det;
        inv[2] = (j[1] * j[5] - j[2] * j[4]) / det;

        inv[3] = (j[5] * j[6] - j[3] * j[8]) / det;
        inv[4] = (j[0] * j[8] - j[2] * j[6]) / det;
        inv[5] = (j[2] * j[3] - j[0] * j[5]) / det;

        inv[6] = (j[3] * j[7] - j[4] * j[6]) / det;
        inv[7] = (j[1] * j[6] - j[0] * j[7]) / det;
        inv[8] = (j[0] * j[4] - j[1] * j[3]) / det;
        return inv;
    }
};

template<int Dim>
class MFSFElement
{
    public:
        explicit MFSFElement(const std::array<int, Dim> &degrees)
            : degrees(degrees), n_basis(ComputeBasisCount(degrees))
        {
            static_assert(Dim == 2 || Dim == 3, "MFSFElement only supports Dim = 2 or 3");
        }

        int GetNumLocalBasis() const
        {
            return n_basis;
        }

        MFSFEvaluation<Dim> EvaluateSingleQP(
            const std::array<std::vector<double>, Dim> &basis_1d,
            const std::array<std::vector<double>, Dim> &dbasis_1d,
            const std::vector<double> &control_points) const
        {
            ValidateInput(basis_1d, dbasis_1d, control_points);

            MFSFEvaluation<Dim> result;
            result.n_basis = n_basis;
            result.basis_values.assign(n_basis, 0.0);
            result.basis_gradients_param.assign(Dim * n_basis, 0.0);
            result.basis_gradients_physical.assign(Dim * n_basis, 0.0);
            result.physical_point.fill(0.0);
            result.jacobian.fill(0.0);
            result.inv_jacobian.fill(0.0);
            result.det_jacobian = 0.0;

            std::array<int, Dim> tensor_index{};
            double partition = 0.0;
            std::array<double, Dim> dpartition{};
            dpartition.fill(0.0);

            for (int basis_id = 0; basis_id < n_basis; ++basis_id)
            {
                DecodeBasisIndex(basis_id, tensor_index);

                double n_value = 1.0;
                for (int axis = 0; axis < Dim; ++axis)
                {
                    n_value *= basis_1d[axis][tensor_index[axis]];
                }

                result.basis_values[basis_id] = n_value;
                partition += n_value;

                for (int deriv_axis = 0; deriv_axis < Dim; ++deriv_axis)
                {
                    double dn_value = 1.0;
                    for (int axis = 0; axis < Dim; ++axis)
                    {
                        if (axis == deriv_axis)
                        {
                            dn_value *= dbasis_1d[axis][tensor_index[axis]];
                        }
                        else
                        {
                            dn_value *= basis_1d[axis][tensor_index[axis]];
                        }
                    }

                    result.basis_gradients_param[deriv_axis * n_basis + basis_id] = dn_value;
                    dpartition[deriv_axis] += dn_value;
                }
            }

            for (int basis_id = 0; basis_id < n_basis; ++basis_id)
            {
                result.basis_values[basis_id] /= partition;
            }

            for (int deriv_axis = 0; deriv_axis < Dim; ++deriv_axis)
            {
                for (int basis_id = 0; basis_id < n_basis; ++basis_id)
                {
                    const double raw_derivative =
                        result.basis_gradients_param[deriv_axis * n_basis + basis_id];
                    result.basis_gradients_param[deriv_axis * n_basis + basis_id] =
                        (raw_derivative - dpartition[deriv_axis] * result.basis_values[basis_id]) / partition;
                }
            }

            for (int basis_id = 0; basis_id < n_basis; ++basis_id)
            {
                for (int coord = 0; coord < Dim; ++coord)
                {
                    const double cp = control_points[basis_id * Dim + coord];
                    result.physical_point[coord] += cp * result.basis_values[basis_id];

                    for (int param_axis = 0; param_axis < Dim; ++param_axis)
                    {
                        result.jacobian[coord * Dim + param_axis] +=
                            cp * result.basis_gradients_param[param_axis * n_basis + basis_id];
                    }
                }
            }

            result.det_jacobian = JacobianOps<Dim>::Determinant(result.jacobian);
            if (Abs(result.det_jacobian) < 1.0e-14)
            {
                throw std::runtime_error("Singular Jacobian in MFSFElement");
            }

            result.inv_jacobian = JacobianOps<Dim>::Inverse(result.jacobian, result.det_jacobian);

            for (int phys_axis = 0; phys_axis < Dim; ++phys_axis)
            {
                for (int basis_id = 0; basis_id < n_basis; ++basis_id)
                {
                    double value = 0.0;
                    for (int param_axis = 0; param_axis < Dim; ++param_axis)
                    {
                        value += result.inv_jacobian[param_axis * Dim + phys_axis]
                            * result.basis_gradients_param[param_axis * n_basis + basis_id];
                    }
                    result.basis_gradients_physical[phys_axis * n_basis + basis_id] = value;
                }
            }

            return result;
        }

    private:
        std::array<int, Dim> degrees;
        int n_basis;

        static int ComputeBasisCount(const std::array<int, Dim> &degrees)
        {
            int count = 1;
            for (int axis = 0; axis < Dim; ++axis)
            {
                count *= (degrees[axis] + 1);
            }
            return count;
        }

        static double Abs(const double value)
        {
            return value >= 0.0 ? value : -value;
        }

        void ValidateInput(const std::array<std::vector<double>, Dim> &basis_1d,
            const std::array<std::vector<double>, Dim> &dbasis_1d,
            const std::vector<double> &control_points) const
        {
            for (int axis = 0; axis < Dim; ++axis)
            {
                const std::size_t expected = static_cast<std::size_t>(degrees[axis] + 1);
                if (basis_1d[axis].size() != expected)
                {
                    throw std::invalid_argument("basis_1d has incorrect axis length");
                }
                if (dbasis_1d[axis].size() != expected)
                {
                    throw std::invalid_argument("dbasis_1d has incorrect axis length");
                }
            }

            if (control_points.size() != static_cast<std::size_t>(n_basis * Dim))
            {
                throw std::invalid_argument("control_points has incorrect size");
            }
        }

        void DecodeBasisIndex(const int flat_index, std::array<int, Dim> &tensor_index) const
        {
            int remainder = flat_index;
            for (int axis = 0; axis < Dim; ++axis)
            {
                const int width = degrees[axis] + 1;
                tensor_index[axis] = remainder % width;
                remainder /= width;
            }
        }
};

#endif
