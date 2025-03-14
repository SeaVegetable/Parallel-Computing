#include "IENGenerator.hpp"

std::vector<int> IENGenerator::GenerateIEN1D(const BSplineBasis * const &basis)
{
    const int p = basis->GetDegree();
    const std::vector<double> S = basis->GetKnotVector();
    const int nFunc = static_cast<int>(S.size()) - p - 1;
    const int nElem = nFunc - p;
    const int nLocBas = p + 1;

    std::vector<int> IEN{};

    int A = 0;
    for (int i = 0; i < nFunc; ++i)
    {
        A++;
        if (i >= nLocBas)
        {
            for (int iloc = p; iloc >= 0; --iloc)
            {
                int B = A - iloc - 1;
                IEN.push_back(B);
            }
        }
    }

    return IEN;
}

std::vector<int> IENGenerator::GenerateIEN2D(
    const BSplineBasis * const &basis1, const BSplineBasis * const &basis2)
{
    const int p = basis1->GetDegree();
    const int q = basis2->GetDegree();
    const std::vector<double> S = basis1->GetKnotVector();
    const std::vector<double> T = basis2->GetKnotVector();
    const int m = static_cast<int>(S.size()) - p - 1;
    const int n = static_cast<int>(T.size()) - q - 1;
    const int nElem = (m - p) * (n - q);
    const int nLocBas = (p + 1) * (q + 1);

    std::vector<int> IEN{};

    int A = 0;
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            A++;
            if (i >= p && j >= q)
            {
                for (int jloc = q; jloc >= 0; --jloc)
                {
                    for (int iloc = p; iloc >= 0; --iloc)
                    {
                        int B = A - iloc - jloc * m - 1;
                        IEN.push_back(B);
                    }
                }
            }
        }
    }

    return IEN;
}