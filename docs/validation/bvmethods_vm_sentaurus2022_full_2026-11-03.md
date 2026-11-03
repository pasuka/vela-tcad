# Sentaurus U-2022.12 compatibility: production A0 complete BV

The user explicitly authorized independent Vela VM AND/OR complete BV scans using production A0/G2 and the explicit `sentaurus_u2022_12` constants/material profile. Each full solver process receives 43200 seconds including cold Poisson; build, regressions and short cold qualification are outside that budget. The existing original-A1 T4 block remains and this task does not launch A1.

## Comparison being followed up

The earlier "14 of 15 J frames align" result compares original native Sentaurus exports with Vela reassembly at exactly the same frozen native states. Five external targets were 0.0100000000, 0.0235533333, 0.0419226178, 0.0668191213 and 0.1005621824 V, with native iterations 1/2/3 at each (15 frames). Existing raw native Jacobian/RHS/directions and physical states were reused; no new Sentaurus simulation was performed. The old CODATA/SRH-fixed and new compatibility Vela observations use the same native states and canonical physical-coordinate mapping, exact prescribed-variable elimination, and shared-state row/column normalization. Metrics are not solver acceptance gates or autonomous-path comparisons.

Normalized maximum RHS difference decreased from 2.79894655e-5 to 9.61177636e-10. Fourteen regular J frames have relative Frobenius differences 4.50470593e-10 to 5.39359308e-8; all Poisson blocks are about 2.1e-15. The first p01_i01 frame remains 1.53909088e-3 overall J difference, localized to six hole rows at nodes229/791/1308/1546/2084/2587; cause remains open. Full BV convergence/workload under the new profile was not established by those frozen observations. See [constants validation](bvmethods_sentaurus2022_constants_2026-11-03.md).

## Execution contract

A fresh source snapshot, independent Linux build and separate control/output directories are used. Compile-time constants select q/epsilon0/kB/h/m0/DOS pi; compatible Formula-1 Nc/Nv/ni are regenerated from the original par with active OldSlotboom. Every other material field, relative dielectric constant, mesh and doping remain unchanged. Binary `--physical-constants` and material profile must agree.

Both groups use automatic derivatives with protected model prediction review and rho reporting, reset each target, retaining original Full candidate handling. This is production A0/G2, without A1 actual-rho review, A1e or selector. AND uses enforce_joint; OR uses enforce; independent load/KCL/nonfinite/other original gates and report-only conservation remain unchanged.

2600 nodes/5020 Tri3, 300K, JD/OldSlotboom, total-impurity mobility, SRH plus GradQF avalanche, BTBT off, R=1e7 ohm um, potential_density/QF SG, floor_adjacent_edges, 1% protection, floor1e-100, exact elimination/log extrapolation/full_step/N40, F1e-8/RMS1/load1e-6 V, original Quasistationary .01 V/Min10uV/Max2500V/Increment1.41/Decrement2; no extra internal-voltage limiter. Each PARDISO process uses one thread.

## Current status

Independent compatibility build exited0 and required features were confirmed. Ten C++groups, materials Python and environment checks all exited0; both cold qualification windows passed. Both full solver processes actually started Beijing2026-11-03 13:46:11, deadline2026-11-04 01:46:11 each. Both groups remain running.

## Terminal qualification and comparison

Success requires solver0, not timed out, |Id|>=.0013 A/um, every accepted point passing its own AND/OR plus all original independent gates, and internal-voltage interpolation at1e-4 A/um within0.5% of6.673833 V. Merely reaching requested external voltage is insufficient. OR residual-only/update-only/both/zero-update counts and report-only carrier/continuity failures are separate. No BV success will be claimed before terminal qualification.

The comparison baseline is the completed CODATA productionA0: AND421/581/160 (accept/attempt/fail), BV6.67030242054 V; OR398/551/153, BV6.67028529897 V. Compare common voltage regions and inclusive attempt+record+other costs. Do not infer sole SRH causation, uniform workload reduction or dedicated-machine speedup from shared-VM wall times.

## Full scan confirmed, 2026-11-03

Independent build exited0 in754.157 s; all10 C++ groups, materials Python and runtime-environment checks exited0. Both cold windows accepted2/2 attempts with0failures and passed their original applicable gates. OR cold branches: residual-only1/both1/update-only0/zero-update1; report carrier/continuity unsatisfied1 each retained.

Both full processes actually started 2026-11-03 13:46:11 Beijing, each deadline 2026-11-04 01:46:11. The43200-second budgets include cold Poisson and exclude all preparation. Binary/core/source identities are saved in the ignored deployment record. At about540seconds, AND had22 accepted/27attempts/4failed, OR25/32/6; IV can lag diagnostic rows. Both recent trials ended in an absolute linear-defect failure after bootstrap andFull retry; this does not invalidate previously accepted states. No full-BV result is claimed while running.

At about1041s, AND29accepted/39attempts/8failed, lastouter6.144558010277/inner5.754732729437V/I3.898254363e-8A/um; OR35/46/10, last4.347782871916/4.221109339988V/I1.266736402e-8A/um. Latest OR target4.513427459421V reports residual-branch convergence but is not yet in the acceptanceCSV. No extra point is counted. Both full scans continue.


## Monitoring snapshot, 2026-11-03 about14:10 Beijing

Both full scans remain running, elapsed about1421s, deadlineNov4 01:46:11. AND40accepted/52attempts/11failed, lastouter7.240641323732/inner6.341850540256V/I8.987908755e-8A/um. OR43/54/11, last5.920342953215/5.583646937590V/I3.366960426e-8A/um. The latest AND7.478846996561V iteration12 row reports F2.255e-10/RMS0.101036/load-1.81288e-8V and convergence, but is not yet in the acceptedCSV; no point is added early. Latest OR accepted under its residual branch, F9.339e-9/RMS560.004/load2.69327e-8V. The terminal Full tag/alpha0/linearNaN does not describe all earlier applied steps or a linear failure.


## Monitoring snapshot, 2026-11-03 about14:31 Beijing

Both full scans remain running after2683s, unchanged deadlineNov4 01:46:11. AND70accepted/93attempts/22failed, lastouter25.47735532665/inner6.659082106905V/I1.881827310e-6A/um. OR73/97/23, last12.65639821172/6.641014553276V/I6.015383659e-7A/um. Latest target28.32951360732V(iter13) and14.55912480081V(iter20) failed the scaled absolute linear-defect requirement after bootstrap+same-stateFull retry:3.615906e-8/2.190359e-8 versus1e-8; componentwise errors passed. These are failed trial states; terminaltrial/alpha0 do not describe all previous steps. NoBV qualification yet.
