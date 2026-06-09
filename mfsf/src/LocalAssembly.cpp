#include "LocalAssembly.hpp"

void LocalAssembly::AssemLocalStiffnessLoad(const Element * const &elem,
    const std::vector<double> &eCP)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> N{};
    std::vector<double> dN_dx{};
    std::vector<double> dN_dy{};
    std::vector<double> invJac{};
    double jacobian;
    double xx, yy;

    ResetStiffnessLoad();

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateElementSingleQP(qp1[ii], qp2[jj], eCP,
                N, dN_dx, dN_dy, xx, yy, jacobian, invJac);

            jacobian_qp.push_back(jacobian);
            invJac_qp.insert(invJac_qp.end(), invJac.begin(), invJac.end());

            double J_W = w1[ii]*w2[jj]*jacobian;

            for (int kk = 0; kk < n; ++kk)
            {
                for (int ll = 0; ll < n; ++ll)
                {
                    Kloc[kk*n+ll] -= J_W * (dN_dx[kk]*dN_dx[ll] + dN_dy[kk]*dN_dy[ll]);
                }
                Floc[kk] += J_W * N[kk] * Getf(xx, yy);
            }
        }
    }
}
