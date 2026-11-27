# Native avalanche current/density identification

The user authorized strict original AvalDensGradQF alignment after the load-line
workflow. That workflow is now qualified through external 37.1 V. This stage
must identify the actual source construction before changing Vela physics.

## Compiler constraint and current experiment

The initial PMI compiler wrapper selected GCC 4.8.5 and failed on its required
`-std=c++14` flag before running the device. The user subsequently supplied an
existing GCC 11 toolchain under `vela_build_20260914/toolchain/opt/rh/devtoolset-11/root/usr`.
It reports GCC 11.2.1 and successfully compiles the diagnostic with unchanged
`cmi -O` options. Sentaurus reports a GCC 9.2.0 build and recommends identical
compiler versions, so compilation alone does not establish runtime compatibility.
Restore the isolated PMI experiment and qualify actual loading and frozen output;
do not install a compiler or alter system libraries. Preserve the initial failure.

The GCC 11 diagnostic now passed actual loading and frozen 10 V output checks
in `build/slot_pmi_gcc11_20260923`: exit 0, unchanged primary state, and exact
nodal alpha/generation parity with the built-in constant-alpha control. The
native launcher supplies its own libstdc++.so.6.0.27; system-shell ldd alone
selects an insufficient older runtime. Coordinate/field matching exposes native
vertex choices, including 17 electron and 9 hole differences from double-precision
maximum-angle selection, all inside the prior near-tie envelope. The vector
magnitude and unmatched traces remain unqualified; see the validation report.

Follow-up now completed: all 15,116 electron and 10,960 hole unique callbacks
map uniquely using coordinate and field checks; prior missing matches were
decimal-rounding keys and over-wide field windows. The AvalDensGradQF-off PMI
control passed exact built-in parity. Full current gradients at electrode-contact
elements resolve the interface-overlap vector discrepancy; non-contact interfaces
retain tangential projection. The inferred current scale and remaining contact
saved-state precision limits are documented in
`docs/validation/slot_ldmos_pmi_callback_followup_2026-09-23.md`.
Correct the isolated research kernel and verify its nine-column derivatives
before coupled integration; no production behavior changed in this follow-up.

The isolated contact correction is now implemented and verified: Release 61
cases/2,702 assertions and 8/8 selected CTests pass. Full nine-column tests retain
live statistics and mobility, including pure-normal contact gradients. Native
10/37 V inputs replay unchanged outside the contact/interface overlap; targeted
contact differences pass the existing 1e-5 derivative-block criterion. See
`docs/validation/slot_ldmos_contact_kernel_fix_2026-09-23.md` for saved-state
precision limits and local source improvements. Coupled assembly remains pending.

The available fallback uses the built-in isotropic van Overstraeten model with
only both carriers' low/high `a` set to 1 cm^-1 and `b` to 0 V/cm. At the fixed
300 K state this should give unit alpha wherever the model evaluates a nonzero
driving field; the output must verify that assumption. Run separate Load/Plot
controls with AvalDensGradQF enabled and disabled, retaining all other physical
parameters. Check state and mobility identity, compare local generation and
integrals, and retain zero-field/contact cases separately. This isolates the
ionization law but does not itself identify an element density or validate a
device solution. Keep the original parameters and a precise parameter diff.

The two constant-alpha controls completed with identical primary states and
mobility/density arrays. Simple harmonic density reconstruction still fails
local source parity. An E0 threshold pair also completed; it revealed a smooth
transition rather than a hard switch. The b=0 coefficient data support a cosine
transition across 0.95--1.05 E0. The completed unit-mobility controls identify
the harmonic mean of minimum/maximum vertex densities in 196 individually
identifiable carrier-elements. A nonzero-b holdout confirms smoothing of final
alpha curves, including the derivative of the mixing weight.

Exact-state vertex-kernel replay identifies local best-vertex IALMob/HFS in the
source current. At 10/37 V, unambiguous interior significant nodes agree below
2.2e-12 relative; the remaining interior differences are covered by explicit
near-equal-angle vertex-choice bounds. The unique native tie rule is not yet
reproduced. Do not replace it with a fitted per-element choice.

Next resolve interface/contact field versus source-current construction and
truncated source volumes. The isolated frozen MixAverageBoxMethod control completed;
its explicit source shares recover constant-alpha output at all 12,410 Si nodes.
Geometric clipping differs by vertex permutation in 57 cells; a universal obtuse
swap rule failed, so preserve imported geometry and its provenance.

An independent C++ research kernel now accepts source shares and selected vertex
explicitly, with full nine-potential chain derivatives. It passed Release unit,
protocol, actual-state derivative and 10/37 V source replay checks. It is not wired
into solver configuration or Newton assembly. The next step is to qualify that
integration after resolving or explicitly bounding the remaining geometry-choice
and weak-source differences. The original-physics element-output control ended:
after fixing a missing parameter input, r2 rejected `eAlphaAvalanche/Element` at
Plot. Do not repeat that unsupported field. Generation/Element was not reached
and remains a separate possible diagnostic; retain the failed logs and partial
TDR as failed evidence, never as a checkpoint.
The measurements and limitations are in
`docs/validation/slot_ldmos_native_avalanche_identification_2026-09-23.md`.

## Restored PMI experiment

Use the installed T-2022.03 public avalanche PMI interface to observe
`currentWoMob`, driving field, and evaluation coordinates. The installed manual
pages 1260–1262 specify the input vector and its units; pages 1240–1241 specify
coordinate access. It is a new measurement, not assumed to equal the generic
runtime `eCurrent/Element` or `ReadScalar(element,eDensity)` data.

Load the same original 10 V saved state, perform no Newton updates, and replace
only the avalanche ionization coefficient with a constant 1 cm^-1 diagnostic
PMI. Its derivatives are identically zero. Preserve all other physical switches,
including AvalDensGradQF, ElementVolumeAvalanche and IALMob, and record the exact
deck diff. Use one native thread for deterministic trace output. Trace electron
and hole calls separately; plot and export the unchanged state and source rates.
This constant-alpha diagnostic is not an accepted device solution and is never
fed into continuation or used as a fitted physical model.

Compare with the identical frozen-state control with only AvalDensGradQF
disabled. Verify state identity, model activation, units, callback support and
generation reconstruction before interpreting any inferred density. In
particular, the manual recommends the direction of currentWoMob for model use;
its magnitude must be independently checked against the constant-alpha source.

If the callback exposes only nodal or otherwise ambiguous data, do not infer an
element formula from it. Prepare an isolated interior-element control, or use
additional independent frozen controls, to disambiguate density and source
distribution. No production formula or coefficient is changed on the basis of
terminal-current proximity.

## Acceptance before implementation

- Reconstruct native local and integrated generation, not only terminal current.
- Identify element density, current/field construction and vertex allocation
  separately, including interior, interface and contact cases.
- Validate a fixed formula at multiple biases without coefficient fitting.
- The standalone explicit-geometry kernel may be checked before all native
  geometry rules are inferred; it must not silently choose those rules. Qualify
  solver integration and its complete assembled Jacobian separately, then
  requalify checkpoints under the new physical configuration.

Native inputs/outputs remain under a new dated diagnostic prefix. Keep all old
reference states and logs unchanged, hash uploaded/downloaded files, and retain
failed experiments. Full BVDS and post-breakdown qualification remain pending.
