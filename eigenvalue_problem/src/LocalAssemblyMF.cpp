#include "LocalAssemblyMF.hpp"

void LocalAssemblyMF::AssemLocalLoad(const Element * const &elem,
    const std::vector<double> &eCP)
{
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const std::vector<double> qp1 = quad1->GetQuadraturePoint();
    const std::vector<double> qp2 = quad2->GetQuadraturePoint();
    const std::vector<double> w1 = quad1->GetWeight();
    const std::vector<double> w2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    std::vector<double> R{};
    std::vector<double> dR_dx{};
    std::vector<double> dR_dy{};
    double jacobian;
    double xx, yy;

    ResetStiffnessLoad();

    for (int jj = 0; jj < nqp2; ++jj)
    {
        for (int ii = 0; ii < nqp1; ++ii)
        {
            elem->GenerateElementSingleQP(qp1[ii], qp2[jj], eCP, R, dR_dx, dR_dy, xx, yy, jacobian);

            double J_W = w1[ii]*w2[jj]*jacobian;

            for (int kk = 0; kk < n; ++kk)
            {
                Floc[kk] += J_W * R[kk] * Getf(xx, yy);
            }
        }
    }
}

void LocalAssemblyMF::LocalMatMulMF(ElementMF * const &elem,
    const std::vector<double> &eCP)
{
    const int n = elem->GetNumLocalBasis();

    std::vector<double> R{};
    std::vector<double> dR_dx{};
    std::vector<double> dR_dy{};
    double jacobian;
    double xx, yy;

    elem->GenerateElement(quad1, quad2, eCP);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const int nqp = nqp1 * nqp2;

    std::vector<double> temp_x(nqp, 0.0);
    std::vector<double> temp_y(nqp, 0.0);

    for (int ii = 0; ii < nqp; ++ii)
    {
        for (int jj = 0; jj < n; ++jj)
        {
            temp_x[ii] += elem->get_dR_dx(ii, jj) * Floc_in[jj];
            temp_y[ii] += elem->get_dR_dy(ii, jj) * Floc_in[jj];
        }
    }

    for (int ii = 0; ii < nqp; ++ii)
    {
        temp_x[ii] *= elem->get_JxW(ii);
        temp_y[ii] *= elem->get_JxW(ii);
    }

    for (int ii = 0; ii < n; ++ii)
    {
        for (int jj = 0; jj < nqp; ++jj)
        {
            Floc_out[ii] += elem->get_R(jj, ii) * temp_x[ii];
            Floc_out[ii] += elem->get_R(jj, ii) * temp_y[ii]; 
        }
    }
}