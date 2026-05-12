#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"

namespace vela {

double sgElectronFlux(double n0, double n1, double dpsi, double Vt,
                      double mu, double h)
{
    // u = (psi_j - psi_i) / Vt
    const double u = dpsi / Vt;
    // J_ij = mu_n * Vt / h * [ B(+u)*n_j - B(-u)*n_i ]
    return mu * Vt / h * (bernoulli(u) * n1 - bernoulli(-u) * n0);
}

double sgHoleFlux(double p0, double p1, double dpsi, double Vt,
                  double mu, double h)
{
    const double u = dpsi / Vt;
    // J_ij = mu_p * Vt / h * [ B(-u)*p_j - B(+u)*p_i ]
    return mu * Vt / h * (bernoulli(-u) * p1 - bernoulli(u) * p0);
}

} // namespace vela
