# Native avalanche current/density identification

The user authorized strict original AvalDensGradQF alignment after the load-line
workflow. That workflow is now qualified through external 37.1 V. This stage
must identify the actual source construction before changing Vela physics.

## Compiler constraint and current experiment

The installed PMI compiler wrapper selects GCC 4.8.5 and fails on its required
`-std=c++14` flag before running the device. Sentaurus reports a GCC 9.2.0 build;
the user confirmed on 2026-09-23 that no compatible compiler is available.
PMI compilation is therefore shelved. Do not install a compiler or weaken ABI
checks as part of this experiment.

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
and weak-source differences. An original-physics element-output control is running
to check whether the native source can expose the selected vertex directly.
The measurements and limitations are in
`docs/validation/slot_ldmos_native_avalanche_identification_2026-09-23.md`.

## Shelved PMI experiment

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
