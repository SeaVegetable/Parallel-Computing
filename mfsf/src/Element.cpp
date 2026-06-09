#include "Element.hpp"

void Element::GenerateElementSingleQP(const double &xi, const double &eta,
    const std::vector<double> &eCP,
    std::vector<double> &N,
    std::vector<double> &dN_dx,
    std::vector<double> &dN_dy,
    double &x, double &y,
    double &jacobian,
    std::vector<double> &invJac) const
{
    BernsteinBasis * bern1 = new BernsteinBasis(p);
    BernsteinBasis * bern2 = new BernsteinBasis(q);

    RefElement * ref = new RefElement();

    const std::vector<double> N1 = ref->GenerateBasis1DSingleQP(bern1, extraction1, xi);
    const std::vector<double> dN1 = ref->GenerateBasisDerivative1DSingleQP(bern1, extraction1, xi, hx);
    const std::vector<double> N2 = ref->GenerateBasis1DSingleQP(bern2, extraction2, eta);
    const std::vector<double> dN2 = ref->GenerateBasisDerivative1DSingleQP(bern2, extraction2, eta, hy);

    delete bern1;
    delete bern2;
    delete ref;

    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    N.clear();

    for (int j = 0; j < q + 1; ++j)
    {
        for (int i = 0; i < p + 1; ++i)
        {
            N.push_back(N1[i] * N2[j]);
            dN_dxi.push_back(dN1[i] * N2[j]);
            dN_deta.push_back(N1[i] * dN2[j]);
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;

    x = 0.0;
    y = 0.0;

    const int n_loc_bas = (p + 1) * (q + 1);
    for (int idx = 0; idx < n_loc_bas; ++idx)
    {
        x += eCP[2 * idx] * N[idx];
        y += eCP[2 * idx + 1] * N[idx];
        dx_dxi += eCP[2 * idx] * dN_dxi[idx];
        dx_deta += eCP[2 * idx] * dN_deta[idx];
        dy_dxi += eCP[2 * idx + 1] * dN_dxi[idx];
        dy_deta += eCP[2 * idx + 1] * dN_deta[idx];
    }

    jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi;

    const double dxi_dx = dy_deta / jacobian;
    const double dxi_dy = -dx_deta / jacobian;
    const double deta_dx = -dy_dxi / jacobian;
    const double deta_dy = dx_dxi / jacobian;

    invJac.clear();
    invJac.push_back(dxi_dx);
    invJac.push_back(deta_dx);
    invJac.push_back(dxi_dy);
    invJac.push_back(deta_dy);

    dN_dx.clear();
    dN_dy.clear();

    for (int idx = 0; idx < n_loc_bas; ++idx)
    {
        dN_dx.push_back(dxi_dx * dN_dxi[idx] + deta_dx * dN_deta[idx]);
        dN_dy.push_back(dxi_dy * dN_dxi[idx] + deta_dy * dN_deta[idx]);
    }

    jacobian *= hx * hy;
}
