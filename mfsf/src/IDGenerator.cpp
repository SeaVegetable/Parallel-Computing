#include "IDGenerator.hpp"

std::vector<int> IDGenerator::GenerateID1D(const int &nFuncX)
{
    std::vector<int> ID(nFuncX, -1);

    for (int i = 1; i < nFuncX-1; ++i) ID[i] = i - 1;  

    return ID;
}

std::vector<int> IDGenerator::GenerateID1D(const BSplineBasis * const &basis)
{
    const int p = basis->GetDegree();
    const std::vector<double> S = basis->GetKnotVector();
    const int nFunc = static_cast<int>(S.size()) - p - 1;

    std::vector<int> ID(nFunc, -1);

    for (int i = 1; i < nFunc-1; ++i) ID[i] = i - 1;  

    return ID;
}

std::vector<int> IDGenerator::GenerateID2D(const int &nFuncX, const int &nFuncY)
{
    std::vector<int> ID(nFuncX*nFuncY, -1);

    int A = nFuncX + 1;
    for (int j = 1; j < nFuncY-1; ++j)
    {
        for (int i = 1; i < nFuncX-1; ++i)
        {
            ID[j*nFuncX+i] = A;
            A++;
        }
        A++;
    }

    return ID;
}

std::vector<int> IDGenerator::GenerateID2D(
    const BSplineBasis * const &basis1, const BSplineBasis * const &basis2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;

    std::vector<int> ID(m*n, -1);

    int A = m + 1;
    for (int j = 1; j < n-1; ++j)
    {
        for (int i = 1; i < m-1; ++i)
        {
            ID[j*m+i] = A;
            A++;
        }
        A++;
    }

    return ID;
}

std::vector<int> IDGenerator::GenerateID3D(const int &nFuncX, const int &nFuncY, const int &nFuncZ)
{
    std::vector<int> ID(nFuncX*nFuncY*nFuncZ, -1);

    int A = nFuncX*nFuncY + nFuncX + 1;
    for (int k = 1; k < nFuncZ-1; ++k)
    {
        for (int j = 1; j < nFuncY-1; ++j)
        {
            for (int i = 1; i < nFuncX-1; ++i)
            {
                ID[k*nFuncX*nFuncY+j*nFuncX+i] = k*nFuncX*nFuncY+j*nFuncX+i;
            }
        }
    }

    return ID;
}

std::vector<int> IDGenerator::GenerateID3D(
    const BSplineBasis * const &basis1,
    const BSplineBasis * const &basis2,
    const BSplineBasis * const &basis3)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const int r = basis3->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const std::vector<double> U = basis3->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;
    const int o = static_cast<int>(U.size()) - r - 1;

    std::vector<int> ID(m*n*o, -1);

    int A = m*n + m + 1;
    for (int k = 1; k < o-1; ++k)
    {
        for (int j = 1; j < n-1; ++j)
        {
            for (int i = 1; i < m-1; ++i)
            {
                ID[k*m*n + j*m + i] = k*m*n + j*m + i;
            }
        }
    }

    return ID;
}