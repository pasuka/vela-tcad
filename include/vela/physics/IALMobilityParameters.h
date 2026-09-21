#pragma once

#include "vela/core/Types.h"

namespace vela {

/** IALMob coefficients in the explicit cm-based units of the source model.
 * These coefficients do not change with the solver's internal unit system.
 * T-2022.03 Device User Guide, equations 290--311 and tables 62--63.
 */
struct IALMobilityCarrierParameters {
    Real mumax = 1417.0; // cm^2/(V s)
    Real mumin = 52.2;
    Real theta = 2.285;
    Real n_ref = 9.68e16; // cm^-3
    Real alpha = 0.68;
    Real ndop_ref = 1.0e18;
    Real nsc_ref = 1.0e18;
    Real S = 0.3042; // K (cm/V)^(2/3)
    Real p = 4.0;
    Real l_crit_c = 1.0e3; // cm
    Real B = 9.0e5; // cm/s
    Real C = 4400.0; // cm^(5/3)/(V^(2/3) s)
    Real lambda = 0.057;
    Real k = 1.0;
    Real alpha_ph2d_A = 1.0;
    Real alpha_ph2d_D = 1.0;
    Real lambda_ph2d_A = 1.0;
    Real lambda_ph2d_D = 1.0;
    Real delta = 3.97e13; // cm^2/(V s), field normalized by 1 V/cm
    Real lambda_sr = 0.057;
    Real A = 2.0;
    Real alpha_sr = 0.0; // cm^3
    Real nu = 0.0;
    Real eta = 1.0e50; // V^2/(cm s)
    Real N1 = 1.0; // cm^-3
    Real N2 = 1.0;
    Real alpha_sr_A = 1.0;
    Real alpha_sr_D = 1.0;
    Real lambda_sr_A = 1.0;
    Real lambda_sr_D = 1.0;
    Real l_crit = 1.0e3; // cm
    Real D1_inv = 135.0; // cm^2/(V s)
    Real D2_inv = 40.0;
    Real nu0_inv = 1.5;
    Real nu1_inv = 2.0;
    Real nu2_inv = 0.5;
    Real alpha1_inv = 0.0;
    Real alpha2_inv = 0.0;
    Real D1_acc = 135.0;
    Real D2_acc = 40.0;
    Real nu0_acc = 1.5;
    Real nu1_acc = 2.0;
    Real nu2_acc = 0.5;
    Real alpha1_acc = 0.0;
    Real alpha2_acc = 0.0;
    Real EnormMinimum = 0.0; // V/cm
};

struct IALMobilityParameters {
    IALMobilityCarrierParameters electron{};
    IALMobilityCarrierParameters hole = [] {
        IALMobilityCarrierParameters value;
        value.mumax = 470.5;
        value.mumin = 44.9;
        value.theta = 2.247;
        value.n_ref = 2.23e17;
        value.alpha = 0.719;
        return value;
    }();
    Real nref_D = 4.0e20; // cm^-3
    Real nref_A = 7.2e20;
    Real cref_D = 0.21;
    Real cref_A = 0.5;
    Real me_over_m0 = 1.0;
    Real mh_over_m0 = 1.258;
    bool full_phumob = true;
    bool clustering_everywhere = false;
    int phonon_combination = 1;
};

} // namespace vela
