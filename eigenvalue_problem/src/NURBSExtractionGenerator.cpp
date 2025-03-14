#include "NURBSExtractionGenerator.hpp"

std::vector<int> NURBSExtractionGenerator::GenerateExtraction1D(const BSplineBasis * const &basis)
{
    const std::vector<double> S = basis->GetKnotVector();
    const int p = basis->GetDegree();
    const int m = static_cast<int>(S.size());
    const int nb = m - p - 1;

    int a = p + 1;
    std::vector<double> C(a * a * nb, 0.0);

    std::vector<double> I(a * a, 0.0);
    for (int i = 0; i < a; ++i) I[i * a + i] = 1.0;

    for (int i = 0; i < a; ++i) C[idxC(i, i, 0, a)] = 1.0;

    int span = 1;
    int A = p + 1;
    int B = A + 1;

    while (B < m)
    {
        for (int i = 0; i < a; ++i) C[idxC(i, i, span+1, a)] = 1.0;
        int iVal = B;
        while (B < m && S[B] == S[iVal])
        {
            B++;
        }
        int mult = B - iVal + 1;
        if (mult < p)
        {
            double numer = S[B - 1] - S[A - 1];

            int r = p - mult;
            std::vector<double> alphas(r, 0.0);

            for (int j = p; j >= mult + 1; --j)
                alphas[j - mult - 1] = numer / (S[A - 1 + j] - S[A - 1]);

            for (int j = 1; j <= r; ++j)
            {
                int save = r - j + 1;
                int s = mult + j;

                for (int k = p + 1; k >= s + 1; --k)
                {
                    double alpha = alphas[k - s - 1];

                    for (int row = 1; row <= a; ++row)
                    {
                        double oldVal_k = C[idxC(row - 1, k - 1, span, a)];
                        double oldVal_km1 = C[idxC(row - 1, k - 2, span, a)];
                        C[idxC(row - 1, k - 1, span + 1, a)] = alpha * oldVal_k + (1.0 - alpha) * oldVal_km1;
                    }
                }

                if (B < m)
                {
                    for (int offset = 0; offset <= j; ++offset)
                    {
                        int rowLeft = save + offset;
                        int rowRight = p - j + offset;

                        C[idxC(rowLeft, rowRight, span + 1, a)] = C[idxC[rowRoght, p + 1, span, a]];
                    }
                }
            }

            span++;
            if (B < m)
            {
                A = B;
                B++;
            }
        }
    }
}