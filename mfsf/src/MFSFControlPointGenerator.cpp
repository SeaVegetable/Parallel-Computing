#include "MFSFControlPointGenerator.hpp"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>

std::vector<double> MFSFControlPointGenerator::SolveLinearSystem(std::vector<double> A,
    std::vector<double> b,
    const int n) const
{
    for (int pivot = 0; pivot < n; ++pivot)
    {
        int max_row = pivot;
        double max_value = A[pivot * n + pivot];
        if (max_value < 0.0) max_value = -max_value;

        for (int row = pivot + 1; row < n; ++row)
        {
            double value = A[row * n + pivot];
            if (value < 0.0) value = -value;
            if (value > max_value)
            {
                max_value = value;
                max_row = row;
            }
        }

        if (max_value < 1.0e-14)
        {
            throw std::runtime_error("Singular matrix in control point generation");
        }

        if (max_row != pivot)
        {
            for (int col = 0; col < n; ++col)
            {
                std::swap(A[pivot * n + col], A[max_row * n + col]);
            }
            std::swap(b[pivot], b[max_row]);
        }

        const double pivot_value = A[pivot * n + pivot];
        for (int row = pivot + 1; row < n; ++row)
        {
            const double factor = A[row * n + pivot] / pivot_value;
            for (int col = pivot; col < n; ++col)
            {
                A[row * n + col] -= factor * A[pivot * n + col];
            }
            b[row] -= factor * b[pivot];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int row = n - 1; row >= 0; --row)
    {
        double sum = b[row];
        for (int col = row + 1; col < n; ++col)
        {
            sum -= A[row * n + col] * x[col];
        }
        x[row] = sum / A[row * n + row];
    }
    return x;
}

std::vector<double> MFSFControlPointGenerator::GenerateControlPoints1D(const BSplineBasis * const &basis,
    const double &min_value,
    const double &max_value) const
{
    if (basis == nullptr)
    {
        throw std::invalid_argument("BSpline basis pointer cannot be null");
    }

    const int p = basis->GetDegree();
    const std::vector<double> S = basis->GetKnotVector();
    const int nFunc = static_cast<int>(S.size()) - p - 1;

    if (nFunc < 2 * p + 2)
    {
        std::vector<double> cp(nFunc, min_value);
        if (nFunc == 1)
        {
            cp[0] = min_value;
            return cp;
        }
        const double step = (max_value - min_value) / static_cast<double>(nFunc - 1);
        for (int i = 0; i < nFunc; ++i)
        {
            cp[i] = min_value + step * i;
        }
        return cp;
    }

    const double ds = (S[p + 1] - S[p]) / (p + 1);
    std::vector<double> sample_points;
    for (int i = 0; i < p; ++i)
    {
        sample_points.push_back(S[p] + (i + 1) * ds);
    }

    const double J = (max_value - min_value) / (S[nFunc] - S[p]);

    std::vector<double> A(p * p, 0.0);
    std::vector<double> b(p, J);
    for (int i = 0; i < p; ++i)
    {
        const int span = basis->FindSpan(sample_points[i]);
        const std::vector<double> dN = basis->DerBasisFuns(sample_points[i], span, 1);
        for (int col = 0; col < p; ++col)
        {
            A[i * p + col] = dN[col + 1];
        }
        b[i] -= min_value * dN[0];
    }

    const std::vector<double> solution = SolveLinearSystem(A, b, p);

    std::vector<double> CP;
    CP.push_back(min_value);
    for (int i = 0; i < p; ++i)
    {
        CP.push_back(solution[i]);
    }

    const int L = nFunc - 2 * p - 2;
    const double d = (max_value + min_value - 2.0 * CP.back()) / (L + 1);
    for (int i = 0; i < L; ++i)
    {
        CP.push_back(CP.back() + d);
    }
    for (int i = p; i > 0; --i)
    {
        CP.push_back(max_value - (CP[i] - min_value));
    }
    CP.push_back(max_value);

    return CP;
}
