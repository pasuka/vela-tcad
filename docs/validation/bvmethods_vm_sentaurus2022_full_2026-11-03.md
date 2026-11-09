# Sentaurus U-2022.12 compatibility: production A0 complete BV

The user explicitly authorized independent Vela VM AND/OR complete BV scans using production A0/G2 and the explicit `sentaurus_u2022_12` constants/material profile. Each full solver process receives 43200 seconds including cold Poisson; build, regressions and short cold qualification are outside that budget. The existing original-A1 T4 block remains and this task does not launch A1.

## Terminal result: both complete BV scans qualified, 2026-11-03

Both production A0/G2 compatibility scans finished before their original 12-hour deadlines. Every accepted point was independently checked against its own original AND/OR and independent gates, both final currents reached 0.0013 A/um, and the 1e-4 A/um internal-voltage interpolation passed the 0.5% reference gate. Reaching the requested 100000 V was not required after the qualified stop current.

| Group | End, Beijing | Solve time h:m:s | Accept/attempt/fail | Final current A/um | BV inner V | Signed reference difference | Interpolation bracket mV |
|---|---|---|---:|---:|---:|---:|---:|
| AND | 2026-11-03 18:06:03.152 | 4:19:51.270 | 395/547/152 | 0.00133815058922 | 6.670290808717 | -0.05307582% | 4.68623707 |
| OR | 2026-11-03 18:03:53.558 | 4:17:41.675 | 407/564/157 | 0.00136389252584 | 6.670284552248 | -0.05316956% | 4.38723485 |

AND report-only carrier/continuity unsatisfied points: 27/8; OR: 60/10. These remain visible and are not newly promoted to enforced gates. OR branches: residual-only143, update-only2, both262; zero-update1 (included in branch counts). All original enforced gates passed.

### Workload comparison against repaired CODATA production A0

| Outer range V | CODATA AND | Compatibility AND | CODATA OR | Compatibility OR |
|---|---:|---:|---:|---:|
| 0-1 | 16/18/2 | 13/13/0 | 21/27/6 | 15/17/2 |
| 1-5 | 15/19/4 | 12/16/4 | 14/16/2 | 23/31/8 |
| 5-10 | 43/59/16 | 25/35/10 | 22/30/8 | 28/38/10 |
| 10-100 | 74/103/29 | 72/100/28 | 71/99/28 | 71/98/27 |
| 100-712.5 | 241/348/107 | 242/350/108 | 240/348/108 | 239/347/108 |
| 712.5-1000 | 16/16/0 | 16/16/0 | 15/15/0 | 16/16/0 |
| 100-1000 | 257/364/107 | 258/366/108 | 255/363/108 | 255/363/108 |
| 1000-3000 | 6/6/0 | 6/6/0 | 7/7/0 | 6/6/0 |

Cells are accept/attempt/direct-fail counts, including the zero-bias anchor in the first interval. The 100-1000 row overlaps the two preceding high-voltage subranges and must not be added to them.

Overall AND changed from421/581/160 to395/547/152 (-26/-34/-8); OR changed from398/551/153 to407/564/157 (+9/+13/+4). Thus there is no uniform workload reduction. Compatibility AND 0-1 V has13/13/0, but this alone does not establish Sentaurus-equivalent stopping semantics or identical targets/states. Native sweep reference is T-2022.03-SP2; the frozen-state/constants target is U-2022.12.

### Common external range and inclusive costs

AND, common outer<= 13388.8589208 V: CODATA 420/580/160, compatibility 395/547/152. Full wall time 19198.989397 -> 15591.270144 s is a shared-VM observation, not a dedicated-machine speedup.

AND compatibility partition: attempt 14847.771534479 + record 738.018593183 + other 5.380079955 = total 15591.170207617 s. CODATA: 18318.303963189 + 863.927855183 + 16.584160660 = 19198.815979032 s. Inclusive substages must not be summed; per-attempt costs were not emitted.

OR, common outer<= 13646.2918922 V: CODATA 397/550/153, compatibility 407/564/157. Full wall time 17344.735478 -> 15461.674915 s is a shared-VM observation, not a dedicated-machine speedup.

OR compatibility partition: attempt 14704.931855614 + record 751.282856954 + other 5.348476461 = total 15461.563189029 s. CODATA: 16501.938603183 + 828.151779598 + 14.502574039 = 17344.592956820 s. Inclusive substages must not be summed; per-attempt costs were not emitted.

### Terminal evidence and identity

Exact supervisor, master worker and both group workers/solvers exited0. No owned task processes remain. All6155 source files and final binary/core/material/input identities were checked; no runtime configuration or deadline was changed. Necessary evidence was downloaded once: ZIP SHA256 `5ed8f806f7de0e397c2c11e1319fdfa96dc8d9747d06049fc19a48988544d126`, 3311373 bytes, 37 data members, all member SHA256 values verified. Successful logs were not read. Private job IDs/paths and raw evidence remain in the ignored task directory. No new simulation was launched for this analysis.

Both scans and the already-qualified T470p window task are finished; the joint monitoring can stop. Original A1 T4 remains blocked. No production stopping-rule change, floor change or further complete BV is authorized by this result.

## Comparison being followed up

The earlier "14 of 15 J frames align" result compares original native Sentaurus exports with Vela reassembly at exactly the same frozen native states. Five external targets were 0.0100000000, 0.0235533333, 0.0419226178, 0.0668191213 and 0.1005621824 V, with native iterations 1/2/3 at each (15 frames). Existing raw native Jacobian/RHS/directions and physical states were reused; no new Sentaurus simulation was performed. The old CODATA/SRH-fixed and new compatibility Vela observations use the same native states and canonical physical-coordinate mapping, exact prescribed-variable elimination, and shared-state row/column normalization. Metrics are not solver acceptance gates or autonomous-path comparisons.

Normalized maximum RHS difference decreased from 2.79894655e-5 to 9.61177636e-10. Fourteen regular J frames have relative Frobenius differences 4.50470593e-10 to 5.39359308e-8; all Poisson blocks are about 2.1e-15. The first p01_i01 frame remains 1.53909088e-3 overall J difference, localized to six hole rows at nodes229/791/1308/1546/2084/2587; cause remains open. Full BV convergence/workload under the new profile was not established by those frozen observations. See [constants validation](bvmethods_sentaurus2022_constants_2026-11-03.md).

## Execution contract

A fresh source snapshot, independent Linux build and separate control/output directories are used. Compile-time constants select q/epsilon0/kB/h/m0/DOS pi; compatible Formula-1 Nc/Nv/ni are regenerated from the original par with active OldSlotboom. Every other material field, relative dielectric constant, mesh and doping remain unchanged. Binary `--physical-constants` and material profile must agree.

Both groups use automatic derivatives with protected model prediction review and rho reporting, reset each target, retaining original Full candidate handling. This is production A0/G2, without A1 actual-rho review, A1e or selector. AND uses enforce_joint; OR uses enforce; independent load/KCL/nonfinite/other original gates and report-only conservation remain unchanged.

2600 nodes/5020 Tri3, 300K, JD/OldSlotboom, total-impurity mobility, SRH plus GradQF avalanche, BTBT off, R=1e7 ohm um, potential_density/QF SG, floor_adjacent_edges, 1% protection, floor1e-100, exact elimination/log extrapolation/full_step/N40, F1e-8/RMS1/load1e-6 V, original Quasistationary .01 V/Min10uV/Max2500V/Increment1.41/Decrement2; no extra internal-voltage limiter. Each PARDISO process uses one thread.

## Current status

Independent compatibility build exited0 and required features were confirmed. Ten C++groups, materials Python and environment checks all exited0; both cold qualification windows passed. Both full solver processes actually started Beijing2026-11-03 13:46:11, deadline2026-11-04 01:46:11 each. Both groups subsequently completed and passed terminal qualification; the terminal result below supersedes the dated running snapshots.

## Terminal qualification and comparison

Success requires solver0, not timed out, |Id|>=.0013 A/um, every accepted point passing its own AND/OR plus all original independent gates, and internal-voltage interpolation at1e-4 A/um within0.5% of6.673833 V. Merely reaching requested external voltage is insufficient. OR residual-only/update-only/both/zero-update counts and report-only carrier/continuity failures are separate. Terminal qualification is now complete; see the qualified result below.

The comparison baseline is the completed CODATA productionA0: AND421/581/160 (accept/attempt/fail), BV6.67030242054 V; OR398/551/153, BV6.67028529897 V. Compare common voltage regions and inclusive attempt+record+other costs. Do not infer sole SRH causation, uniform workload reduction or dedicated-machine speedup from shared-VM wall times.

## Full scan confirmed, 2026-11-03

Independent build exited0 in754.157 s; all10 C++ groups, materials Python and runtime-environment checks exited0. Both cold windows accepted2/2 attempts with0failures and passed their original applicable gates. OR cold branches: residual-only1/both1/update-only0/zero-update1; report carrier/continuity unsatisfied1 each retained.

Both full processes actually started 2026-11-03 13:46:11 Beijing, each deadline 2026-11-04 01:46:11. The43200-second budgets include cold Poisson and exclude all preparation. Binary/core/source identities are saved in the ignored deployment record. At about540seconds, AND had22 accepted/27attempts/4failed, OR25/32/6; IV can lag diagnostic rows. Both recent trials ended in an absolute linear-defect failure after bootstrap andFull retry; this does not invalidate previously accepted states. No full-BV result is claimed while running.

At about1041s, AND29accepted/39attempts/8failed, lastouter6.144558010277/inner5.754732729437V/I3.898254363e-8A/um; OR35/46/10, last4.347782871916/4.221109339988V/I1.266736402e-8A/um. Latest OR target4.513427459421V reports residual-branch convergence but is not yet in the acceptanceCSV. No extra point is counted. Both full scans continue.


## Monitoring snapshot, 2026-11-03 about14:10 Beijing

Both full scans remain running, elapsed about1421s, deadlineNov4 01:46:11. AND40accepted/52attempts/11failed, lastouter7.240641323732/inner6.341850540256V/I8.987908755e-8A/um. OR43/54/11, last5.920342953215/5.583646937590V/I3.366960426e-8A/um. The latest AND7.478846996561V iteration12 row reports F2.255e-10/RMS0.101036/load-1.81288e-8V and convergence, but is not yet in the acceptedCSV; no point is added early. Latest OR accepted under its residual branch, F9.339e-9/RMS560.004/load2.69327e-8V. The terminal Full tag/alpha0/linearNaN does not describe all earlier applied steps or a linear failure.


## Monitoring snapshot, 2026-11-03 about14:31 Beijing

Both full scans remain running after2683s, unchanged deadlineNov4 01:46:11. AND70accepted/93attempts/22failed, lastouter25.47735532665/inner6.659082106905V/I1.881827310e-6A/um. OR73/97/23, last12.65639821172/6.641014553276V/I6.015383659e-7A/um. Latest target28.32951360732V(iter13) and14.55912480081V(iter20) failed the scaled absolute linear-defect requirement after bootstrap+same-stateFull retry:3.615906e-8/2.190359e-8 versus1e-8; componentwise errors passed. These are failed trial states; terminaltrial/alpha0 do not describe all previous steps. NoBV qualification yet.

## Monitoring snapshot, 2026-11-03 about14:46 Beijing

Both full scans remain running after3606s, deadlineNov4 01:46:11 unchanged. AND90accepted/121attempts/30failed, lastouter42.082643097705/inner6.656823509707 V/I3.54258195979e-6 A/um. OR97/132/33, lastouter34.017705391719/inner6.658287008128 V/I2.73594184865e-6 A/um. LatestAND42.996654390739 V(iter6,F1.84416e-10/RMS.039739/load1.20189e-7) andOR35.603305376921 V(retry1,iter5,F3.30995e-9/RMS12.9888/load-6.83458e-8) report convergence but are not yet acceptedCSV rows. No point is counted early; terminalalpha0/linearNaN does not classify earlier actual updates. NoBV result is claimed. The separate T470p window pair has ended and been qualified, leaving only these VM full scans under joint monitoring.

## Monitoring snapshot, 2026-11-03 about15:01 Beijing

Both full scans remain running after4488s; deadlineNov4 01:46:11 unchanged. AND108accepted/147attempts/37failed, lastouter70.314427678936/inner6.652116092371 V/I6.366231166646e-6 A/um. OR120/164/42, lastouter63.616546983612/inner6.653081888972 V/I5.696346505007e-6 A/um. LatestAND73.860243039066 V(iter5,F3.4465e-10/RMS.188687/load1.78949e-7) reports convergence but is not yet in acceptedCSV. LatestOR71.464676116851 V(iter5) failed the scaled absolute linear defect after bootstrap and same-stateFull retry:5.47223e-4 versus1e-8; componentwise9.04067e-17 passed1e-10. This is a failed trial state, not an accepted state. Terminaltrial/alpha0/zero protection do not describe all earlier applied steps or an actual zero update. NoBV success is claimed; T470p remains completed and is not queried.

## Monitoring snapshot, 2026-11-03 15:16 Beijing

AND: full_scan, elapsed5389.9s; accepted/attempts/direct-failures=129/177/46; last accepted outer114.62023853034221/inner6.6476884286582587 V/I1.0797255016549812e-05 A/um. Latest diagnostic target120.19132177821324 V/retry1/iteration5, F2.4370420463547806e-10/RMS0.19877633932515226/load-6.8449651280388935e-08 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

OR: full_scan, elapsed5389.9s; accepted/attempts/direct-failures=141/194/51; last accepted outer110.04241895720796/inner6.6480219009915134 V/I1.033943966846121e-05 A/um. Latest diagnostic target118.27228671917005 V/retry0/iteration40, F0.0073067398754498693/RMS237051.27177742642/load-0.0018039790804493805 V; terminal label=trial/observing_model_failure, review(model/linear)=0/0, protection(n/p)=0/0, converged0; failure=coupled Newton iteration budget exhausted.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.

## Monitoring snapshot, 2026-11-03 15:31 Beijing

AND: full_scan, elapsed6314.3s; accepted/attempts/direct-failures=149/205/54; last accepted outer169.93577506313727/inner6.6450605571644159 V/I1.6329071464033688e-05 A/um. Latest diagnostic target175.95926917967293 V/retry0/iteration6, F1.8397569285465339e-10/RMS0.00052472294714109898/load8.07562798854633e-08 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

OR: full_scan, elapsed6314.3s; accepted/attempts/direct-failures=164/225/60; last accepted outer171.11470208421585/inner6.6450262606733705 V/I1.6446967599205039e-05 A/um. Latest diagnostic target180.40686360560974 V/retry0/iteration10, F14.902493340139518/RMS2549245947.3050499/load-0.23105286773358102 V; terminal label=trial/bootstrap_linear_failure_restore, review(model/linear)=0/1, protection(n/p)=0/0, converged0; failure=bootstrap linear solve failed: direct bordered refinement did not meet componentwise accuracy: 3.463043e-07; full Jacobian retry failed: direct bordered refined residual exceeds the scaled tolerance.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.

## Monitoring snapshot, 2026-11-03 15:46 Beijing

AND: full_scan, elapsed7198.2s; accepted/attempts/direct-failures=170/235/63; last accepted outer233.09368786853483/inner6.6440520959080152 V/I2.2644963565450508e-05 A/um. Latest diagnostic target239.65722736751346 V/retry0/iteration6, F3.4755865599025027e-10/RMS0.0002135821443217458/load1.3236331142252311e-07 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

OR: full_scan, elapsed7198.2s; accepted/attempts/direct-failures=182/253/68; last accepted outer228.23266865216004/inner6.6440795946397611 V/I2.2158858914440304e-05 A/um. Latest diagnostic target237.52883698712591 V/retry0/iteration5, F4.3567509800139141e-10/RMS0.020528582203294052/load-1.501733208897349e-07 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.

## Monitoring snapshot, 2026-11-03 16:00 Beijing

AND: full_scan, elapsed8080.7s; accepted/attempts/direct-failures=193/267/73; last accepted outer308.99326864917902/inner6.6443662801398711 V/I3.023489022397645e-05 A/um. Latest diagnostic target311.72595341931213 V/retry1/iteration5, F3.1624664185468756e-10/RMS0.003124632275957012/load-1.053148253049585e-07 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

OR: full_scan, elapsed8080.6s; accepted/attempts/direct-failures=206/285/78; last accepted outer307.55851040068211/inner6.6443496143242138 V/I3.0091416105403399e-05 A/um. Latest diagnostic target310.06202335935114 V/retry1/iteration5, F3.9602871515803088e-10/RMS0.00068559414706118916/load1.6342966091542621e-07 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.

## Monitoring snapshot, 2026-11-03 16:15 Beijing

AND: full_scan, elapsed8983.0s; accepted/attempts/direct-failures=213/296/82; last accepted outer378.54897174927333/inner6.6455389952102406 V/I3.7190343268985354e-05 A/um. Latest diagnostic target388.62872553108872 V/retry0/iteration6, F4071.4195673960667/RMS209747719.24396434/load-11921.222372620394 V; terminal label=trial/bootstrap_linear_failure_restore, review(model/linear)=0/1, protection(n/p)=0/0, converged0; failure=bootstrap linear solve failed: direct bordered refined residual exceeds the scaled tolerance; full Jacobian retry failed: direct bordered refined residual exceeds the scaled tolerance.

OR: full_scan, elapsed8982.8s; accepted/attempts/direct-failures=226/317/88; last accepted outer376.89527370515577/inner6.6455038275930862 V/I3.702497699250967e-05 A/um. Latest diagnostic target392.67213421578924 V/retry0/iteration40, F687.36857128783504/RMS69080821854.455765/load-7.3249792451684357 V; terminal label=trial/observing_model_failure, review(model/linear)=0/0, protection(n/p)=0/0, converged0; failure=coupled Newton iteration budget exhausted.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.

## Monitoring snapshot, 2026-11-03 16:30 Beijing

AND: full_scan, elapsed9885.5s; accepted/attempts/direct-failures=236/329/92; last accepted outer455.83573160701241/inner6.6474913734063934 V/I4.4918824023054427e-05 A/um. Latest diagnostic target466.02535667273736 V/retry0/iteration13, F4090.9017903300919/RMS191918796300.57477/load-305.31187141890103 V; terminal label=trial/bootstrap_linear_failure_restore, review(model/linear)=0/1, protection(n/p)=0/0, converged0; failure=bootstrap linear solve failed: direct bordered refined residual exceeds the scaled tolerance; full Jacobian retry failed: direct bordered refined residual exceeds the scaled tolerance.

OR: full_scan, elapsed9885.5s; accepted/attempts/direct-failures=252/352/98; last accepted outer467.21827623430636/inner6.6478234075747906 V/I4.6057045299261002e-05 A/um. Latest diagnostic target473.36752074513805 V/retry0/iteration5, F2.3885192589229219e-10/RMS0.020386391947728694/load-1.9814382312688394e-07 V; terminal label=full/full_step_retained, review(model/linear)=0/0, protection(n/p)=0/0, converged1; failure=none.

Latest diagnostic rows are trial/terminal snapshots, not additional acceptedCSV points. Terminalalpha0/linearNaN/zero protection and mode labels do not establish earlier actual updates or matrix history. Deadline2026-11-04 01:46:11 Beijing is unchanged; noBV terminal qualification yet. T470p is finished and not queried.
