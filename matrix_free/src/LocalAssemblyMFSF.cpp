#include "LocalAssemblyMFSF.hpp"

void LocalAssemblyMFSF::AssemLocalLoad(ElementMF * const &elem,
    const std::vector<double> &eCP)
{
    elem->GenerateElement(quad1, quad2, eCP);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const int nqp = nqp1 * nqp2;
    const int n = elem->GetNumLocalBasis();

    ResetLoad();

    for (int ii = 0; ii < nqp; ++ii)
    {
        double x = 0.0;
        double y = 0.0;
        for (int jj = 0; jj < n; ++jj)
        {
            x += elem->get_R(ii, jj) * eCP[2*jj];
            y += elem->get_R(ii, jj) * eCP[2*jj+1];
        }

        for (int jj = 0; jj < n; ++jj)
        {
            Floc[jj] += elem->get_R(ii, jj) * Getf(x, y) * elem->get_JxW(ii);
        }
    }
}

void LocalAssemblyMFSF::LocalMatMulMF(ElementMF * const &elem,
    const std::vector<double> &eCP)
{
    elem->GenerateElement(quad1, quad2, eCP);
    const int nqp1 = quad1->GetNumQuadraturePoint();
    const int nqp2 = quad2->GetNumQuadraturePoint();
    const double qw1 = quad1->GetWeight();
    const double qw2 = quad2->GetWeight();
    const int n = elem->GetNumLocalBasis();

    ResetStiffnessLoadOut();

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int ii = 0; ii < nqp1; ++ii)
            {
                for (int jj = 0; jj < nqp2; ++jj)
                {

                }
            }
        }
    }
}