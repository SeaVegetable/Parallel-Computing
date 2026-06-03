#ifndef MFSF_PRECOMPUTE_HPP
#define MFSF_PRECOMPUTE_HPP

#include <array>
#include <stdexcept>
#include <vector>

#include "BernsteinBasis.hpp"
#include "BSplineBasis.hpp"
#include "MFSFConnectivityGenerators.hpp"
#include "MFSFControlPointGenerator.hpp"
#include "MFSFElement.hpp"
#include "NURBSExtractionGenerator.hpp"
#include "QuadraturePoint.hpp"
#include "RefElement.hpp"

template<int Dim>
struct MFSFPrecomputeData
{
    std::array<int, Dim> degrees;
    std::array<int, Dim> num_funcs;
    std::array<int, Dim> num_elements;
    std::array<int, Dim> nqp;
    std::vector<int> id;
    std::vector<int> ien;
    std::vector<double> control_points;
    std::array<std::vector<double>, Dim> element_sizes;
    std::array<std::vector<double>, Dim> quadrature_points;
    std::array<std::vector<double>, Dim> quadrature_weights;
    std::vector<double> inv_jacobian;
    std::vector<double> det_jacobian;
};

template<int Dim>
struct MFSFQuadratureArrayBuilder;

template<>
struct MFSFQuadratureArrayBuilder<2>
{
    static std::array<QuadraturePoint, 2> Build(const std::array<int, 2> &degrees)
    {
        return {{
            QuadraturePoint(degrees[0] + 1, 0.0, 1.0),
            QuadraturePoint(degrees[1] + 1, 0.0, 1.0)
        }};
    }
};

template<>
struct MFSFQuadratureArrayBuilder<3>
{
    static std::array<QuadraturePoint, 3> Build(const std::array<int, 3> &degrees)
    {
        return {{
            QuadraturePoint(degrees[0] + 1, 0.0, 1.0),
            QuadraturePoint(degrees[1] + 1, 0.0, 1.0),
            QuadraturePoint(degrees[2] + 1, 0.0, 1.0)
        }};
    }
};

template<int Dim>
struct MFSFBernsteinArrayBuilder;

template<>
struct MFSFBernsteinArrayBuilder<2>
{
    static std::array<BernsteinBasis, 2> Build(const std::array<int, 2> &degrees)
    {
        return {{
            BernsteinBasis(degrees[0]),
            BernsteinBasis(degrees[1])
        }};
    }
};

template<>
struct MFSFBernsteinArrayBuilder<3>
{
    static std::array<BernsteinBasis, 3> Build(const std::array<int, 3> &degrees)
    {
        return {{
            BernsteinBasis(degrees[0]),
            BernsteinBasis(degrees[1]),
            BernsteinBasis(degrees[2])
        }};
    }
};

template<int Dim>
class MFSFPrecompute
{
    public:
        MFSFPrecomputeData<Dim> BuildFromBSpline(
            const std::array<const BSplineBasis *, Dim> &bases,
            const std::array<double, Dim> &min_corner,
            const std::array<double, Dim> &max_corner) const
        {
            std::array<int, Dim> degrees{};
            std::array<int, Dim> num_funcs{};
            std::array<int, Dim> num_elements{};
            std::array<std::vector<double>, Dim> knot_vectors;
            std::array<std::vector<double>, Dim> extraction_all;
            std::array<std::vector<double>, Dim> element_sizes;

            NURBSExtractionGenerator extraction_generator;
            MFSFControlPointGenerator control_point_generator;
            MFSFIENGenerator<Dim> ien_generator;
            MFSFIDGenerator<Dim> id_generator;

            for (int axis = 0; axis < Dim; ++axis)
            {
                if (bases[axis] == nullptr)
                {
                    throw std::invalid_argument("BSpline basis pointer cannot be null");
                }
                degrees[axis] = bases[axis]->GetDegree();
                num_funcs[axis] = bases[axis]->GetNumFunctions();
                num_elements[axis] = bases[axis]->GetNumElements();
                knot_vectors[axis] = bases[axis]->GetKnotVector();
                extraction_all[axis] = extraction_generator.GenerateExtraction1D(bases[axis]);
                element_sizes[axis] = ComputeElementSizes(knot_vectors[axis], degrees[axis], num_elements[axis]);
            }

            MFSFPrecomputeData<Dim> data;
            data.degrees = degrees;
            data.num_funcs = num_funcs;
            data.num_elements = num_elements;
            data.id = id_generator.GenerateFromBSpline(bases);
            data.ien = ien_generator.GenerateFromBSpline(bases);
            data.control_points = control_point_generator.GenerateControlPoints<Dim>(bases, min_corner, max_corner);
            data.element_sizes = element_sizes;

            std::array<QuadraturePoint, Dim> quadratures = MFSFQuadratureArrayBuilder<Dim>::Build(degrees);
            for (int axis = 0; axis < Dim; ++axis)
            {
                data.nqp[axis] = quadratures[axis].GetNumQuadraturePoint();
                data.quadrature_points[axis] = quadratures[axis].GetQuadraturePoint();
                data.quadrature_weights[axis] = quadratures[axis].GetWeight();
            }

            const int total_elements = Product(num_elements);
            const int total_qp = Product(data.nqp);
            const int n_loc_bas = Product(LocalWidths(degrees));

            data.inv_jacobian.assign(total_elements * total_qp * Dim * Dim, 0.0);
            data.det_jacobian.assign(total_elements * total_qp, 0.0);

            MFSFElement<Dim> element(degrees);
            RefElement ref;
            std::array<BernsteinBasis, Dim> bernstein = MFSFBernsteinArrayBuilder<Dim>::Build(degrees);

            std::array<int, Dim> elem_index{};
            std::array<int, Dim> qp_index{};
            std::array<int, Dim> local_index{};
            std::array<int, Dim> local_width = LocalWidths(degrees);

            const int extraction_offsets[1] = {0};
            (void)extraction_offsets;

            for (int elem_id = 0; elem_id < total_elements; ++elem_id)
            {
                DecodeIndex(elem_id, num_elements, elem_index);

                std::array<std::vector<double>, Dim> local_extraction{};
                for (int axis = 0; axis < Dim; ++axis)
                {
                    const int block_size = local_width[axis] * local_width[axis];
                    const int offset = elem_index[axis] * block_size;
                    local_extraction[axis].assign(extraction_all[axis].begin() + offset,
                        extraction_all[axis].begin() + offset + block_size);
                }

                std::vector<double> local_control_points(n_loc_bas * Dim, 0.0);
                for (int local_basis_id = 0; local_basis_id < n_loc_bas; ++local_basis_id)
                {
                    DecodeIndex(local_basis_id, local_width, local_index);
                    std::array<int, Dim> global_index{};
                    for (int axis = 0; axis < Dim; ++axis)
                    {
                        global_index[axis] = elem_index[axis] + local_index[axis];
                    }
                    const int global_basis_id = FlattenIndex(global_index, num_funcs);
                    for (int coord = 0; coord < Dim; ++coord)
                    {
                        local_control_points[local_basis_id * Dim + coord] =
                            data.control_points[global_basis_id * Dim + coord];
                    }
                }

                for (int qp_id = 0; qp_id < total_qp; ++qp_id)
                {
                    DecodeIndex(qp_id, data.nqp, qp_index);

                    std::array<std::vector<double>, Dim> basis_1d{};
                    std::array<std::vector<double>, Dim> dbasis_1d{};

                    for (int axis = 0; axis < Dim; ++axis)
                    {
                        const double xi = data.quadrature_points[axis][qp_index[axis]];
                        const std::vector<double> B = bernstein[axis].GetBernsteinBasisSingleQP(xi);
                        const std::vector<double> dB = bernstein[axis].GetBernsteinBasisDerivativeSingleQP(xi);
                        basis_1d[axis] = ref.GenerateBasis1DSingleQP(B, local_extraction[axis]);
                        dbasis_1d[axis] = ref.GenerateBasisDerivative1DSingleQP(dB, local_extraction[axis]);
                        for (std::size_t i = 0; i < dbasis_1d[axis].size(); ++i)
                        {
                            dbasis_1d[axis][i] /= data.element_sizes[axis][elem_index[axis]];
                        }
                    }

                    const MFSFEvaluation<Dim> eval =
                        element.EvaluateSingleQP(basis_1d, dbasis_1d, local_control_points);

                    data.det_jacobian[elem_id * total_qp + qp_id] = eval.det_jacobian;
                    const int jac_offset = (elem_id * total_qp + qp_id) * Dim * Dim;
                    for (int i = 0; i < Dim * Dim; ++i)
                    {
                        data.inv_jacobian[jac_offset + i] = eval.inv_jacobian[i];
                    }
                }
            }

            return data;
        }

    private:
        static std::array<int, Dim> LocalWidths(const std::array<int, Dim> &degrees)
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

        static std::vector<double> ComputeElementSizes(const std::vector<double> &knots,
            const int degree,
            const int num_elements)
        {
            std::vector<double> sizes;
            sizes.reserve(num_elements);
            for (int elem = 0; elem < num_elements; ++elem)
            {
                sizes.push_back(knots[degree + elem + 1] - knots[degree + elem]);
            }
            return sizes;
        }

};

#endif
