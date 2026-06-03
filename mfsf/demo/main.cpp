#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include "BSplineBasis.hpp"
#include "MFSFConnectivityGenerators.hpp"
#include "MFSFElement.hpp"
#include "MFSFH5Writer.hpp"
#include "MFSFPrecompute.hpp"

template<int Dim>
void PrintEvaluation(const MFSFEvaluation<Dim> &eval)
{
    std::cout << "physical_point:";
    for (int i = 0; i < Dim; ++i)
    {
        std::cout << " " << eval.physical_point[i];
    }
    std::cout << std::endl;

    std::cout << "jacobian:";
    for (int i = 0; i < Dim * Dim; ++i)
    {
        std::cout << " " << eval.jacobian[i];
    }
    std::cout << std::endl;

    std::cout << "inv_jacobian:";
    for (int i = 0; i < Dim * Dim; ++i)
    {
        std::cout << " " << eval.inv_jacobian[i];
    }
    std::cout << std::endl;

    std::cout << "det_jacobian: " << eval.det_jacobian << std::endl;
}

void Run2D()
{
    MFSFElement<2> element({1, 1});

    std::array<std::vector<double>, 2> basis_1d = {{
        {0.5, 0.5},
        {0.5, 0.5}
    }};
    std::array<std::vector<double>, 2> dbasis_1d = {{
        {-1.0, 1.0},
        {-1.0, 1.0}
    }};
    std::vector<double> control_points = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        1.0, 1.0
    };

    const MFSFEvaluation<2> eval = element.EvaluateSingleQP(basis_1d, dbasis_1d, control_points);
    std::cout << "[2D]" << std::endl;
    PrintEvaluation(eval);
}

void Run3D()
{
    MFSFElement<3> element({1, 1, 1});

    std::array<std::vector<double>, 3> basis_1d = {{
        {0.5, 0.5},
        {0.5, 0.5},
        {0.5, 0.5}
    }};
    std::array<std::vector<double>, 3> dbasis_1d = {{
        {-1.0, 1.0},
        {-1.0, 1.0},
        {-1.0, 1.0}
    }};
    std::vector<double> control_points = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0,
        1.0, 1.0, 1.0
    };

    const MFSFEvaluation<3> eval = element.EvaluateSingleQP(basis_1d, dbasis_1d, control_points);
    std::cout << "[3D]" << std::endl;
    PrintEvaluation(eval);
}

void RunConnectivityExamples()
{
    std::cout << "[Connectivity 2D]" << std::endl;
    BSplineBasis bx(1, {0.0, 0.0, 0.5, 1.0, 1.0});
    BSplineBasis by(1, {0.0, 0.0, 0.5, 1.0, 1.0});

    MFSFIENGenerator<2> ien2d;
    MFSFIDGenerator<2> id2d;
    const std::vector<int> ien_values_2d = ien2d.GenerateFromBSpline({&bx, &by});
    const std::vector<int> id_values_2d = id2d.GenerateFromBSpline({&bx, &by});

    std::cout << "IEN size: " << ien_values_2d.size() << std::endl;
    std::cout << "ID size: " << id_values_2d.size() << std::endl;
    std::cout << "First 16 IEN:";
    for (std::size_t i = 0; i < ien_values_2d.size() && i < 16; ++i)
    {
        std::cout << " " << ien_values_2d[i];
    }
    std::cout << std::endl;
    std::cout << "ID:";
    for (std::size_t i = 0; i < id_values_2d.size(); ++i)
    {
        std::cout << " " << id_values_2d[i];
    }
    std::cout << std::endl;

    std::cout << "[Connectivity 3D]" << std::endl;
    BSplineBasis bz(1, {0.0, 0.0, 1.0, 1.0});
    MFSFIENGenerator<3> ien3d;
    MFSFIDGenerator<3> id3d;
    const std::vector<int> ien_values_3d = ien3d.GenerateFromBSpline({&bx, &by, &bz});
    const std::vector<int> id_values_3d = id3d.GenerateFromBSpline({&bx, &by, &bz});

    std::cout << "IEN size: " << ien_values_3d.size() << std::endl;
    std::cout << "ID size: " << id_values_3d.size() << std::endl;
    std::cout << "First 24 IEN:";
    for (std::size_t i = 0; i < ien_values_3d.size() && i < 24; ++i)
    {
        std::cout << " " << ien_values_3d[i];
    }
    std::cout << std::endl;
    std::cout << "ID:";
    for (std::size_t i = 0; i < id_values_3d.size(); ++i)
    {
        std::cout << " " << id_values_3d[i];
    }
    std::cout << std::endl;
}

void RunPrecompute2D()
{
    std::cout << "[Precompute 2D]" << std::endl;
    BSplineBasis bx(1, {0.0, 0.0, 0.5, 1.0, 1.0});
    BSplineBasis by(1, {0.0, 0.0, 0.5, 1.0, 1.0});

    MFSFPrecompute<2> precompute;
    const MFSFPrecomputeData<2> data =
        precompute.BuildFromBSpline({&bx, &by}, {0.0, 0.0}, {1.0, 1.0});
    MFSFH5Writer writer;
    writer.Write("mfsf_precompute_2d.h5", data);

    std::cout << "control_points size: " << data.control_points.size() << std::endl;
    std::cout << "ID size: " << data.id.size() << ", IEN size: " << data.ien.size() << std::endl;
    std::cout << "det_jacobian size: " << data.det_jacobian.size() << std::endl;
    std::cout << "inv_jacobian size: " << data.inv_jacobian.size() << std::endl;
    std::cout << "wrote: mfsf_precompute_2d.h5" << std::endl;
    std::cout << "first detJ:";
    for (std::size_t i = 0; i < data.det_jacobian.size() && i < 8; ++i)
    {
        std::cout << " " << data.det_jacobian[i];
    }
    std::cout << std::endl;
    std::cout << "first 8 invJ entries:";
    for (std::size_t i = 0; i < data.inv_jacobian.size() && i < 8; ++i)
    {
        std::cout << " " << data.inv_jacobian[i];
    }
    std::cout << std::endl;
}

void RunPrecompute3D()
{
    std::cout << "[Precompute 3D]" << std::endl;
    BSplineBasis bx(1, {0.0, 0.0, 0.5, 1.0, 1.0});
    BSplineBasis by(1, {0.0, 0.0, 0.5, 1.0, 1.0});
    BSplineBasis bz(1, {0.0, 0.0, 0.5, 1.0, 1.0});

    MFSFPrecompute<3> precompute;
    const MFSFPrecomputeData<3> data =
        precompute.BuildFromBSpline({&bx, &by, &bz}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    MFSFH5Writer writer;
    writer.Write("mfsf_precompute_3d.h5", data);

    std::cout << "control_points size: " << data.control_points.size() << std::endl;
    std::cout << "ID size: " << data.id.size() << ", IEN size: " << data.ien.size() << std::endl;
    std::cout << "det_jacobian size: " << data.det_jacobian.size() << std::endl;
    std::cout << "inv_jacobian size: " << data.inv_jacobian.size() << std::endl;
    std::cout << "wrote: mfsf_precompute_3d.h5" << std::endl;
    std::cout << "first detJ:";
    for (std::size_t i = 0; i < data.det_jacobian.size() && i < 8; ++i)
    {
        std::cout << " " << data.det_jacobian[i];
    }
    std::cout << std::endl;
    std::cout << "first 9 invJ entries:";
    for (std::size_t i = 0; i < data.inv_jacobian.size() && i < 9; ++i)
    {
        std::cout << " " << data.inv_jacobian[i];
    }
    std::cout << std::endl;
}

int main()
{
    std::cout << std::fixed << std::setprecision(6);
    Run2D();
    Run3D();
    RunConnectivityExamples();
    RunPrecompute2D();
    RunPrecompute3D();
    return 0;
}
