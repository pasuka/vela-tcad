# PN2D model and input audit

This campaign compares the actual shipped parameter laws and exported fields.
Matching model names alone is not evidence. The native parameter file is
`../source/models.par`, with SHA256
`b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742`.
The independent T-2022.03-SP2 parameter dump and active-section comparison are
recorded in [parameter_audit.json](results/parameter_audit.json).

## Mobility

At 300 K, P0/P2/P3/D1/D2/G1/G2 use Masetti Formula 1:

`mu(N) = mu_min1 exp(-Pc/N) + (mu_const-mu_min2)/(1+(N/Cr)^alpha) - mu1/(1+(Cs/N)^beta)`.

| Parameter | Electron | Hole |
|---|---:|---:|
| mu_const [cm2/V/s] | 1417 | 470.5 |
| mu_min1 [cm2/V/s] | 52.2 | 44.9 |
| mu_min2 [cm2/V/s] | 52.2 | 0 |
| mu1 [cm2/V/s] | 43.4 | 29.0 |
| Pc [cm^-3] | 0 | 9.23e16 |
| Cr [cm^-3] | 9.68e16 | 2.23e17 |
| Cs [cm^-3] | 3.43e20 | 6.10e20 |
| alpha | 0.68 | 0.719 |
| beta | 2 | 2 |

Vela uses the existing `cell_reconstructed_total_impurity` evaluation policy.
Its law implementation is `src/physics/MobilityModel.cpp`; field comparisons
check the resulting node mobilities, including the compensated junction plane.
P1 disables mobility dependence in SDevice and selects constant mobility in Vela;
both explicitly use 1417/470.5 cm2/V/s. High-field and surface terms are disabled.

## SRH and intrinsic density

SRH uses `R=(np-ni_eff^2)/(tau_p(n+ni_eff)+tau_n(p+ni_eff))`, with a midgap trap,
tau_n=1e-5 s and tau_p=3e-6 s. Neither deck selects doping-dependent lifetime.
Vela's implementation is `src/physics/RecombinationModel.cpp`. P2 disables the
mechanism in both inputs; the exported source integral must also vanish within
the frozen absolute floor.

OldSlotboom uses Ebgn=0.009 eV, Nref=1e17 cm^-3, C=0.5 and
`DeltaEg=Ebgn[ln(N/Nref)+sqrt(ln(N/Nref)^2+C)]`.
The active DD path evaluates N=ND+NA, and
`ni_eff=ni exp(DeltaEg/(2 Vt))`. Specifically,
`include/vela/equation/AssemblerUtils.h::buildEffectiveNodeNi` passes zero
carrier arguments to the generic law in `src/physics/BandgapNarrowing.cpp`.
The generic class can evaluate max(total impurity,n,p), but that does not mean
the active transport ni is updated from injected carrier populations in this
Boltzmann campaign. Fermi-statistics correction is off.
Effective intrinsic density and SRH fields are directly compared at the five
saved biases, so agreement is not inferred from this formula description.

The generator evaluates the shipped silicon bandgap and DOS laws at 300 K.
It includes the OldSlotboom bandgap reference correction dEg0=-0.01595 eV in
the base ni. P3 removes both the narrowing law and that associated reference
correction. This is required to disable the complete model convention; keeping
the corrected ni after disabling BGN would leave part of the model active.
No intrinsic density or lifetime is adjusted to match a terminal current.

## Contacts, doping and units

Both contacts are ideal Ohmic contacts. In Boltzmann equilibrium, charge
neutrality and mass action require n-p=ND-NA and np=ni_eff^2. Vela computes
the majority root using `hypot` and obtains the minority population by division
to avoid subtraction cancellation; see `src/physics/CarrierStatistics.cpp`.
The contact quasi-Fermi potentials equal the applied terminal voltage. The
electrostatic offset follows the same local intrinsic-density convention.
Full spatial comparisons include the contact nodes.
The additional [contact-state audit](results/contact_state_audit.json) also
reports every electrode node, including minority populations excluded by the
global carrier mask. In the completed stage-A states, mass action and neutrality
hold to approximately 1e-14 relative and contact quasi-Fermi errors are below
4e-16 V. The largest minority-carrier difference is about 0.0544%, consistent
with the roughly 0.0272% effective-ni difference through the same mass-action
relation. These raw diagnostic values do not introduce a new threshold or fit
the intrinsic-density constants.

The source SDE uses inclusive adjacent doping windows. At x=1 um, the exported
node has both donor and acceptor concentrations; the campaign preserves this
reported compensation in both solvers. D1/D2 export fresh nodal concentrations,
including the changed junction-plane values. Region-average doping cannot
override this file. For each stage-A mesh, P1/P2/P3 copy the exact P0 binary TDR.

SDevice uses AreaFactor=1 and 1 um out-of-plane width. Port comparisons are in
A/um. Vela's raw electron and hole terminal columns are q times particle inflow;
conventional total current is electron minus hole, so the raw hole column is
negated for comparison to native hCurrent. A conservative cut toward +x equals
minus the Cathode terminal current. No empirical width or sign factor is fitted.

Both branches solve equilibrium independently in Vela, then continue Vela's own
accepted states. Frozen own-state diagnostics recover fields and SG face fluxes;
their exit codes and currents are never used as independent-solve evidence.

## Zero-source continuity diagnostics

The solver's `global_*_continuity_closure_ratio` divides the mismatch by the
largest NET contact flux, integrated source and configured source floor. The
current schema qualifies that ratio only when the integrated source reaches
the source floor. With P2's exactly zero SRH source, roundoff divided by
roundoff can therefore equal one while physical through-current conservation
remains accurate. Applying a 1% cutoff blindly to this column is incorrect.

The comparator verifies finite numeric ratios, enforces the 1% ratio for
source-qualified points, and independently checks each carrier's net terminal
current at unqualified points against the frozen 1% through-current tolerance
and 1e-22 A/um absolute floor. The two-terminal total KCL gate remains 0.1%.
At the five saved states, both tools' port components are additionally checked
against their spatially integrated SRH source. No solver tolerance or contract
threshold is changed to handle the zero-source case.

## Local-contact current-vector recovery evidence

At E1/0.8 V, the legacy nodal conductivity-times-quasi-Fermi-gradient recovery
has total-vector weighted L2 errors of 61.24% (G1) and 85.09% (G2).
Applying that identical recovery to the frozen exported Sentaurus state gives
61.27% and 85.12%. Thus this difference is primarily an output recovery issue;
it is not evidence of a correspondingly wrong independently solved current.
No state, model, default, port current, or acceptance threshold was changed.
Existing dual-face SG recovery gives 0.509%/0.756% errors; existing cell-first
area-weighted SG recovery gives 0.323%/0.450%. The latter has about 2.48% error
in the prescribed contact-endpoint regions and below 0.19% outside them.
Its Jy/Jx norms are 0.11526/0.19010, against native 0.11523/0.19006.
The original field and errors remain available; SG reconstructions are additional
frozen-state evidence and still do not expose Sentaurus internal face fluxes.
