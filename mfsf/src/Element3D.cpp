#include "Element3D.hpp"

void Element3D::GenerateElementSingleQP(
    const double &xi, const double &eta, const double &zeta,
    const std::vector<double> &eCP,
    std::vector<double> &N,
    std::vector<double> &dN_dx,
    std::vector<double> &dN_dy,
    std::vector<double> &dN_dz,
    double &x, double &y, double &z,
    double &jacobian,
    std::vector<double> &invJac) const
{
    BernsteinBasis * bern1 = new BernsteinBasis(p);
    BernsteinBasis * bern2 = new BernsteinBasis(q);
    BernsteinBasis * bern3 = new BernsteinBasis(r);

    RefElement * ref = new RefElement();

    const std::vector<double> N1 = ref->GenerateBasis1DSingleQP(bern1, extraction1, xi);
    const std::vector<double> dN1 = ref->GenerateBasisDerivative1DSingleQP(bern1, extraction1, xi, hx);
    const std::vector<double> N2 = ref->GenerateBasis1DSingleQP(bern2, extraction2, eta);
    const std::vector<double> dN2 = ref->GenerateBasisDerivative1DSingleQP(bern2, extraction2, eta, hy);
    const std::vector<double> N3 = ref->GenerateBasis1DSingleQP(bern3, extraction3, zeta);
    const std::vector<double> dN3 = ref->GenerateBasisDerivative1DSingleQP(bern3, extraction3, zeta, hz);

    delete bern1;
    delete bern2;
    delete bern3;
    delete ref;

    std::vector<double> dN_dxi{};
    std::vector<double> dN_deta{};
    std::vector<double> dN_dzeta{};
    N.clear();

    for (int k = 0; k < r + 1; ++k)
    {
        for (int j = 0; j < q + 1; ++j)
        {
            for (int i = 0; i < p + 1; ++i)
            {
                N.push_back(N1[i] * N2[j] * N3[k]);
                dN_dxi.push_back(dN1[i] * N2[j] * N3[k]);
                dN_deta.push_back(N1[i] * dN2[j] * N3[k]);
                dN_dzeta.push_back(N1[i] * N2[j] * dN3[k]);
            }
        }
    }

    double dx_dxi = 0.0;
    double dx_deta = 0.0;
    double dx_dzeta = 0.0;
    double dy_dxi = 0.0;
    double dy_deta = 0.0;
    double dy_dzeta = 0.0;
    double dz_dxi = 0.0;
    double dz_deta = 0.0;
    double dz_dzeta = 0.0;

    x = 0.0;
    y = 0.0;
    z = 0.0;

    const int n_loc_bas = (p + 1) * (q + 1) * (r + 1);
    for (int idx = 0; idx < n_loc_bas; ++idx)
    {
        x += eCP[3 * idx] * N[idx];
        y += eCP[3 * idx + 1] * N[idx];
        z += eCP[3 * idx + 2] * N[idx];

        dx_dxi += eCP[3 * idx] * dN_dxi[idx];
        dx_deta += eCP[3 * idx] * dN_deta[idx];
        dx_dzeta += eCP[3 * idx] * dN_dzeta[idx];

        dy_dxi += eCP[3 * idx + 1] * dN_dxi[idx];
        dy_deta += eCP[3 * idx + 1] * dN_deta[idx];
        dy_dzeta += eCP[3 * idx + 1] * dN_dzeta[idx];

        dz_dxi += eCP[3 * idx + 2] * dN_dxi[idx];
        dz_deta += eCP[3 * idx + 2] * dN_deta[idx];
        dz_dzeta += eCP[3 * idx + 2] * dN_dzeta[idx];
    }

    jacobian =
        dx_dxi * (dy_deta * dz_dzeta - dy_dzeta * dz_deta) -
        dx_deta * (dy_dxi * dz_dzeta - dy_dzeta * dz_dxi) +
        dx_dzeta * (dy_dxi * dz_deta - dy_deta * dz_dxi);

    const double dxi_dx = (dy_deta * dz_dzeta - dy_dzeta * dz_deta) / jacobian;
    const double dxi_dy = (dx_dzeta * dz_deta - dx_deta * dz_dzeta) / jacobian;
    const double dxi_dz = (dx_deta * dy_dzeta - dx_dzeta * dy_deta) / jacobian;

    const double deta_dx = (dy_dzeta * dz_dxi - dy_dxi * dz_dzeta) / jacobian;
    const double deta_dy = (dx_dxi * dz_dzeta - dx_dzeta * dz_dxi) / jacobian;
    const double deta_dz = (dx_dzeta * dy_dxi - dx_dxi * dy_dzeta) / jacobian;

    const double dzeta_dx = (dy_dxi * dz_deta - dy_deta * dz_dxi) / jacobian;
    const double dzeta_dy = (dx_deta * dz_dxi - dx_dxi * dz_deta) / jacobian;
    const double dzeta_dz = (dx_dxi * dy_deta - dx_deta * dy_dxi) / jacobian;

    invJac.clear();
    invJac.push_back(dxi_dx);
    invJac.push_back(deta_dx);
    invJac.push_back(dzeta_dx);
    invJac.push_back(dxi_dy);
    invJac.push_back(deta_dy);
    invJac.push_back(dzeta_dy);
    invJac.push_back(dxi_dz);
    invJac.push_back(deta_dz);
    invJac.push_back(dzeta_dz);

    dN_dx.clear();
    dN_dy.clear();
    dN_dz.clear();

    for (int idx = 0; idx < n_loc_bas; ++idx)
    {
        dN_dx.push_back(
            dxi_dx * dN_dxi[idx] +
            deta_dx * dN_deta[idx] +
            dzeta_dx * dN_dzeta[idx]);
        dN_dy.push_back(
            dxi_dy * dN_dxi[idx] +
            deta_dy * dN_deta[idx] +
            dzeta_dy * dN_dzeta[idx]);
        dN_dz.push_back(
            dxi_dz * dN_dxi[idx] +
            deta_dz * dN_deta[idx] +
            dzeta_dz * dN_dzeta[idx]);
    }

    jacobian *= hx * hy * hz;
}
