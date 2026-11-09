# Sentaurus U-2022.12 constants compatibility, 2026-11-03

The user authorized a separate compatibility profile. CODATA2018 remains the default; the new profile is used only for newly built validation binaries and regenerated materials. The existing VM T2/T3 experiments retain their frozen cores. No new full BV is started by this task.

## Implementation

`VELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12` selects the audited q, epsilon0, kB, h, m0 and DOS pi in the core. The frontend material emitter takes the same explicit profile and regenerates Formula-1 Nc/Nv and ni while retaining Varshni/OldSlotboom semantics. Relative dielectric constants, solver policies and acceptance gates are unchanged. Built-in no-BGN defaults have a separate reproducible snapshot. The constants profile is immutable per binary.

The material root and compatibility deck declare the expected profile; mismatches fail before assembly. Untagged external materials retain CODATA2018 semantics. Default CODATA deck output stays identical to the existing goldens. `--physical-constants` exposes a binary's actual values. See [configuration](../config_schema.md#physical-constants-profiles) and [audited native values](bvmethods_native_parameter_constants_audit_2026-11-03.md).

At 300 K with the BV Formula-1 parameters and active OldSlotboom:

| Quantity | Compatibility value |
| --- | --- |
| Nc | 2.8582520073866568e19 cm^-3 |
| Nv | 3.1046289324748386e19 cm^-3 |
| base ni | 14638805412.559162 cm^-3 |
| kBT/q | 0.025851995266484913 V |

No fitted Poisson/DOS coefficients are included. DOS pi affects DOS only; geometric/trigonometric pi and hbar keep mathematical pi. The independent nonlocal BTBT source-constant option remains explicit; this BV validation has BTBT off.

## Completed verification

T470p configured and built the independent compatibility profile successfully;
all 6154 snapshot files were checked individually before and after execution.
The binary's `--physical-constants` values match the audited profile. The worker
exited 0 and no owned solver/compiler/worker processes remained. The 318-file
compact evidence package was downloaded once and every member SHA256 checked.
Detailed job paths and receipts remain in the ignored task directory.

| Check | Result |
| --- | --- |
| Material emitter Python | 38 tests passed on local host and T470p |
| Compatibility C++ | constants/material contracts, SRH, carrier statistics, AV contract, coupled load line, PARDISO all passed |
| Default constants / built-in materials | 2 cases / 10 assertions and 3 cases / 192 assertions passed |
| Frozen operators | six packed states finite in F, Full J and bootstrap J |
| Native same-state observations | 15 finite J/RHS frames, 15 raw directions available |
| Native post-state Full continuation | all five passed strict AND and independent gates after 1 update |
| Cold AND, 0 to 0.01 V | 2 accepted / 2 attempts / 0 failures; original gates passed |
| Cold OR, 0 to 0.01 V | 2 accepted / 2 attempts / 0 failures; own OR and independent gates passed |

The cold OR branches were residual-only 1, both 1, update-only 0, with one
zero-update starting point. Each cold group retained one report-only carrier
and one report-only continuity failure; these were not converted into enforced
gates or hidden. All use the existing 2600-node/5020-Tri3 mesh, 300 K,
JD/OldSlotboom, total-impurity mobility, SRH plus GradQF avalanche, BTBT off,
R=1e7 ohm um, and one MKL PARDISO thread. No Sentaurus simulation was run.
The six old packed states are finite-operator observations under the new
constants, not claims of replaying the old nonlinear trajectory.

The default deck regression initially detected seven metadata-only golden
changes. Default deck payload was restored; compatibility decks alone carry
an explicit new profile. The subsequent 75 existing checks had no remaining
product failure (7 unavailable-runner checks skipped); the newly added profile
propagation check passed after correcting its fixture filename. This frontend
metadata correction and its test are the only two files changed after the
T470p snapshot; no C++/physics/material-emitter code changed. The correction
was verified locally. The complete default Release build subsequently exited 0. Full CTest
registered 1219 tests: 1213 passed, 6 skipped, 0 failed (227.446 s); the
controller exited 0 and the JUnit XML and skip list are preserved. This also
covered the corrected frontend profile propagation test. Constants-task
qualification is complete.

## J/RHS comparison

The comparison reuses the original canonical physical-coordinate mapping,
exact prescribed-variable elimination and diagonal/shared-state normalization.
Both the old SRH-fixed observations and new profile use the same 15 native
states. These are normalization-dependent comparison metrics, not new solver
gates or a proof of identical nonlinear paths.

| Metric | Previous CODATA profile | Compatibility profile |
| --- | --- | --- |
| Maximum normalized RHS difference, across 15 frames | 2.79894655e-5 | 9.61177636e-10 |
| Whole-J relative Frobenius, 14 regular frames | previously within 4.2647e-4 to 7.9744e-4 across the full set | 4.50470593e-10 to 5.39359308e-8 |
| Poisson block relative Frobenius, all 15 frames | coefficient bias present | 2.09e-15 to 2.14e-15 |
| First frame p01_i01 whole-J relative Frobenius | see previous per-frame evidence | 1.53909088e-3; retained outlier |
| Final native-iteration frames: one-Full update RMS | 1.26964 to 1.37340 | 1.99871e-4 to 0.308978 |

The true native post-state continuations are a separate experiment: their
single-update RMS values were 1.60589e-5, 1.24852e-4, 5.55066e-4, 1.67636e-3,
and 3.03402e-3; all ended with F below 6e-11 and passed the original strict
AND/load/KCL/other gates. The earlier CODATA post-state experiment required
two updates. An update RMS below one does not establish bit-identical roots.

The remaining first-frame discrepancy is concentrated in six hole rows at
nodes 229, 791, 1308, 1546, 2084 and 2587. Only those six rows exceed 1e-8
in the normalized entry difference; some off-diagonal signs differ. Their
carrier population is about 191 to 266 cm^-3 in n+ regions. The observed AV
source dominates the very small local residual (approximately 1e-25 to
1e-23), while SRH is about 1e-29 to 1e-28. This identifies a localized,
near-zero-current tangent discrepancy worth checking separately; it does not prove
its tangent cause. In particular, no production AV change or smoothing was
made to hide this outlier. The maximum error across all 15 J frames has not
improved, despite the strong improvement in the other 14 frames and RHS.

## Scope and remaining work

The constants mismatch explains the main systematic Poisson/statistical RHS
bias in these observations. It does not remove every J discrepancy or prove
vendor implementation identity. Full BV accuracy and workload under the new
profile have not been established; this task ran short cold windows only.
The existing T4 admission block and running VM T2/T3 cores remain unchanged.
The complete default regression passed as recorded above; the remaining VM
T2/T3 diagnostics continue under their existing authorization.


Two real oneMKL checks skipped in the main default suite were subsequently run
with the installed runtime explicitly selected: both passed, 38 assertions.
The original JUnit counts are retained (1213 passed, 6 skipped); four fixture
E2E checks remain unexecuted in that full-suite invocation.
