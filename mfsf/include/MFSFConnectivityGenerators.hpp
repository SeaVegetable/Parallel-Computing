#ifndef MFSF_CONNECTIVITY_GENERATORS_HPP
#define MFSF_CONNECTIVITY_GENERATORS_HPP

#include <array>
#include <stdexcept>
#include <vector>

#include "BSplineBasis.hpp"

template<int Dim>
class MFSFIENGenerator
{
    public:
        std::vector<int> GenerateFromNumFuncs(const std::array<int, Dim> &num_funcs,
            const std::array<int, Dim> &degrees) const
        {
            Validate(num_funcs, degrees);

            const std::array<int, Dim> num_elem = ComputeNumElem(num_funcs, degrees);
            const std::array<int, Dim> local_width = ComputeLocalWidths(degrees);
            const int n_elem = Product(num_elem);
            const int n_loc_bas = Product(local_width);

            std::vector<int> ien;
            ien.reserve(n_elem * n_loc_bas);

            std::array<int, Dim> elem_index{};
            std::array<int, Dim> local_index{};
            std::array<int, Dim> global_index{};

            for (int elem_id = 0; elem_id < n_elem; ++elem_id)
            {
                DecodeIndex(elem_id, num_elem, elem_index);

                for (int local_basis_id = 0; local_basis_id < n_loc_bas; ++local_basis_id)
                {
                    DecodeIndex(local_basis_id, local_width, local_index);
                    for (int axis = 0; axis < Dim; ++axis)
                    {
                        global_index[axis] = elem_index[axis] + local_index[axis];
                    }
                    ien.push_back(FlattenIndex(global_index, num_funcs));
                }
            }

            return ien;
        }

        std::vector<int> GenerateFromBSpline(const std::array<const BSplineBasis *, Dim> &bases) const
        {
            std::array<int, Dim> num_funcs{};
            std::array<int, Dim> degrees{};
            for (int axis = 0; axis < Dim; ++axis)
            {
                if (bases[axis] == nullptr)
                {
                    throw std::invalid_argument("BSpline basis pointer cannot be null");
                }
                num_funcs[axis] = bases[axis]->GetNumFunctions();
                degrees[axis] = bases[axis]->GetDegree();
            }
            return GenerateFromNumFuncs(num_funcs, degrees);
        }

    private:
        static void Validate(const std::array<int, Dim> &num_funcs,
            const std::array<int, Dim> &degrees)
        {
            for (int axis = 0; axis < Dim; ++axis)
            {
                if (degrees[axis] < 0)
                {
                    throw std::invalid_argument("degree must be non-negative");
                }
                if (num_funcs[axis] <= degrees[axis])
                {
                    throw std::invalid_argument("num_funcs must be larger than degree in every axis");
                }
            }
        }

        static std::array<int, Dim> ComputeNumElem(const std::array<int, Dim> &num_funcs,
            const std::array<int, Dim> &degrees)
        {
            std::array<int, Dim> num_elem{};
            for (int axis = 0; axis < Dim; ++axis)
            {
                num_elem[axis] = num_funcs[axis] - degrees[axis];
            }
            return num_elem;
        }

        static std::array<int, Dim> ComputeLocalWidths(const std::array<int, Dim> &degrees)
        {
            std::array<int, Dim> widths{};
            for (int axis = 0; axis < Dim; ++axis)
            {
                widths[axis] = degrees[axis] + 1;
            }
            return widths;
        }

        static int Product(const std::array<int, Dim> &values)
        {
            int product = 1;
            for (int axis = 0; axis < Dim; ++axis)
            {
                product *= values[axis];
            }
            return product;
        }

        static void DecodeIndex(const int flat,
            const std::array<int, Dim> &widths,
            std::array<int, Dim> &multi_index)
        {
            int remainder = flat;
            for (int axis = 0; axis < Dim; ++axis)
            {
                multi_index[axis] = remainder % widths[axis];
                remainder /= widths[axis];
            }
        }

        static int FlattenIndex(const std::array<int, Dim> &multi_index,
            const std::array<int, Dim> &widths)
        {
            int flat = 0;
            int stride = 1;
            for (int axis = 0; axis < Dim; ++axis)
            {
                flat += multi_index[axis] * stride;
                stride *= widths[axis];
            }
            return flat;
        }
};

template<int Dim>
class MFSFIDGenerator
{
    public:
        std::vector<int> GenerateFromNumFuncs(const std::array<int, Dim> &num_funcs) const
        {
            for (int axis = 0; axis < Dim; ++axis)
            {
                if (num_funcs[axis] < 2)
                {
                    throw std::invalid_argument("num_funcs must be at least 2 in every axis");
                }
            }

            const int total_size = Product(num_funcs);
            std::vector<int> id(total_size, -1);

            std::array<int, Dim> multi_index{};
            int equation = 0;
            for (int flat = 0; flat < total_size; ++flat)
            {
                DecodeIndex(flat, num_funcs, multi_index);

                bool on_boundary = false;
                for (int axis = 0; axis < Dim; ++axis)
                {
                    if (multi_index[axis] == 0 || multi_index[axis] == num_funcs[axis] - 1)
                    {
                        on_boundary = true;
                        break;
                    }
                }

                if (!on_boundary)
                {
                    id[flat] = equation;
                    ++equation;
                }
            }

            return id;
        }

        std::vector<int> GenerateFromBSpline(const std::array<const BSplineBasis *, Dim> &bases) const
        {
            std::array<int, Dim> num_funcs{};
            for (int axis = 0; axis < Dim; ++axis)
            {
                if (bases[axis] == nullptr)
                {
                    throw std::invalid_argument("BSpline basis pointer cannot be null");
                }
                num_funcs[axis] = bases[axis]->GetNumFunctions();
            }
            return GenerateFromNumFuncs(num_funcs);
        }

    private:
        static int Product(const std::array<int, Dim> &values)
        {
            int product = 1;
            for (int axis = 0; axis < Dim; ++axis)
            {
                product *= values[axis];
            }
            return product;
        }

        static void DecodeIndex(const int flat,
            const std::array<int, Dim> &widths,
            std::array<int, Dim> &multi_index)
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
