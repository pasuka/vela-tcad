# VM near-zero-current AV tangents and native stopping audit, 2026-11-03

**Terminal result:** all authorized bounded observations and evidence checks are complete. The observed native update-error branch gives the same0-100 V numerical path as original AND:123 accepted points,165 attempts,42 failures and1448 actual writebacks. It does not reduce work in this controlled scan. This is not a test of every Sentaurus stopping branch: RHSMin remains unqualified, and the near-zero AV native internal convention is not uniquely identified. All production identities and acceptance gates remain unchanged. Earlier progress entries below are historical; the final qualification at the end supersedes them.

The user authorized further near-zero-current AV derivative checks and the proposed sequence of stopping-definition audit, shadow evaluation, controlled same-state targets and bounded adaptive comparison. This task uses an independent diagnostic wrapper linked to the already-qualified `sentaurus_u2022_12` core. Production source, core library, constants/materials, physical models and defaults are unchanged. No new Sentaurus or complete BV scan is launched. The original A1 admission block remains.

## Configuration and evidence boundaries

The device has2600 nodes/5020 Tri3, 300 K, JD/OldSlotboom statistics, total-impurity mobility, SRH plus GradQF avalanche, BTBT off, R=1e7 ohm um, population coordinates/QF SG, floor_adjacent_edges, floor1e-100 and1% population protection. The original A0 automatic/Full handling, N40, independent load/KCL/nonfinite checks and report-only conservation are retained. Each solver uses one PARDISO thread; at most two Newton solver processes run concurrently. Source snapshot has6155 files; private job/path records are in the ignored task directory.

The original native whole-sweep reference is T-2022.03-SP2; the15 frozen states and constants target are U-2022.12. These version scopes are not interchangeable. The present study uses the existing native frozen data, not a new vendor solve.

## AV derivative qualification

Twelve native-state local rows were exported: hole nodes229/791/1308/1546/2084/2587 in p01_i01, plus the same six nodes in p01_i02. All source terms and local production Jacobian rows were retained. The independently reconstructed double transport and AV values closed their production rows before high-precision evaluation. The diagnostic exporter links the unchanged compatibility core.

The250-digit nonlinear kernel uses four positive population perturbation scales,1e-12/1e-20/1e-32/1e-48, and checks every local supporting column. Two explicit arithmetic contracts are compared: (a) keep the exported rounded QF base values, differentiating their physical population dependence locally; (b) reconstruct QF from the same psi/n/p and material statistics in multiprecision. These are diagnostic arithmetic alternatives, not changes to the production model.

| Sample | AD vs250-digit FD with rounded QF anchor | Change after physical high-precision QF reconstruction |
|---|---:|---:|
| First native frame, six hole rows | weighted relative error7.92e-16 to1.16e-15 | about0.2144% to0.8795% |
| Next native iteration, same six rows | weighted relative error5.09e-17 to2.41e-15 | about2.28e-8 to1.49e-7 |

The maximum local QF displacement is approximately1.82e-16 V in the first frame and1.40e-16 V in the next. No local AD/FD sign disagreement was observed in either arithmetic contract; refinement differences were explicitly retained. Thus production AV AD is consistent with its rounded constitutive base values in these rows. Near-zero-current arithmetic sensitivity exists, but that alone does not explain the vendor discrepancy.

### Frozen AV matrix counterfactual

The six VM-measured complete rows replace the corresponding rows of the existing compatibility observation. Exact prescribed-variable elimination, fixed equation-unit conversion from an independent stable Jacobian, one fixed native diagonal scaling and shared physical column scales are used in every arm. Only AV tangents change; source values and states are fixed, and no Newton step is taken.

| Six-row AV tangent arm | Whole-J relative Frobenius difference from native first frame |
|---|---:|
| Production AD |1.53909088096e-3 |
|250-digit derivative at rounded QF anchor |1.53909088096e-3 |
|250-digit derivative after physical QF reconstruction |1.53999449753e-3 |
| AV tangent zero, source value retained |7.72010509021e-4 |
| AV tangent sign reflected, source value retained |1.60359172020e-5 |

The observed first-frame discrepancy is largely opposite to the current AV tangent. This is a frozen diagnostic direction, not authorization or justification to reverse a production derivative. The current AD already differentiates its own residual correctly; remaining uncertainty includes the native near-zero current/norm convention, arithmetic orientation and unexported intermediate values. Next-iteration signed-tangent controls have completed. With the same fixed normalization, its production whole-J difference is1.34787699e-7, zero-AV tangent worsens it to6.80866153e-4 and sign-reflected AV tangent to1.36173127e-3. This sharply limits the sign observation to the first near-zero-current frame; a global sign reversal would damage the next-state agreement. The additional10-file next-state evidence packet is1422317 bytes, SHA256 `36952f63832ea8b0f13c08ff5e430840ae40c13c3336da1e4de42d35205a6d80`, all member hashes verified. Its row-unit/diagonal convention differs from the earlier per-operator-diagonal historical comparison, so its unmodified1.35e-7 number must not be silently substituted into that historical metric. No global AV sign change or smoothing was made.

The initial offline counterfactual omitted equation-row unit conversion; it is preserved as unqualified R1. Corrected R2 reproduces the established first-frame1.53909e-3 baseline before any new arm is interpreted. No MP or Newton execution was repeated for that postprocessing correction.

## Stopping definitions and qualifications

Installed U-2022.12 manual pages204-205 describe the default residual-or-update rule and explicit AND/RhsMinFactor/CheckRhsAfterUpdate alternatives. The current deck does not explicitly request these overrides. Prior native observations establish factor-one raw-direction numerator, post-protection denominator, weighted RMS, Digits5, carrier ErrRef1e10 cm^-3, and metric counts2600+2121+2121+11. These counts include11 aliases of the contact voltage and differ from the number of independent linear unknowns.

The actual current compatibility C++ metric was evaluated using15 original native pre-states, original native raw directions and actual native post-states. All15 results fall inside their native log printing intervals, both with the current rounded potential/contact reference0.025852 V and an exact-VT reference. The maximum relative difference between these two reference choices is4.22e-8. The first target's three reconstructed current-reference errors are2750.96677691,80.2377827407 and0.249964176497; the native values print as2.75e3,8.02e1 and2.4996e-1. This compares the same native directions rather than confusing direction differences with metric differences.

### RHS branch limitation

The native exported RHS infinity norms reproduce the printed iteration norms. Conversion of the circuit column requires physical V, distinct from the packed normalized potential coordinate. A fixed per-equation row-unit conversion identified from one stable native/Vela Jacobian was held against the other frozen states without fitting RHS values or window outcomes. Correcting the contact-column conversion brings the maximum norm relative discrepancy to about2.47e-5 across these frames; however, remaining absolute vector differences are larger than the native1e-5 threshold near roots. These observations do not qualify a complete native RHS branch. They remain shadow-only, and Vela's F threshold is not simply set to1e-5.

The actual experimental third arm is explicitly `native_update`: the verified observed native update-error branch, retaining independent physical gates. It excludes an unqualified native RHS escape branch and therefore is not claimed to implement every Sentaurus stopping semantic. Original AND and OR remain separate controls. First-entry/no-history behavior and the zero anchor remain explicit.

## Controlled execution and current results

The new wrapper builds independently. With intervention disabled, both cold AND and OR pairs pass their original gates and match original-entry IV/updates/HDF datasets bit for bit. Four isolated numerical contract checks pass, covering disabled parity, update-only stopping despite a tighter residual criterion, no false pass without a valid accepted update, and the zero-anchor case. Test R1 used an invalid zero-iteration fixture, failed before a meaningful numerical check, and was retained; R2 uses a valid budget and an explicit failed-linear-solve no-writeback negative control. No runtime iteration budget was changed to repair the fixture.

Eighteen fixed-target controls have completed: six frozen starting states at the first five native external targets, each under original AND, original OR and native_update. Within each triple the initial packed state and initial RHS hashes match. This selected set has identical stopping iteration counts in the three arms. It covers early native targets through0.100562 V; it is not a sample of all difficult0.5-1 V attempts.

| Cold adaptive0-1 V arm | Accept/attempt/direct fail | Applicable original independent gates | Actual solve time s |
|---|---:|---|---:|
| Original AND |14/14/0 | passed |208.197 |
| Original OR |16/18/2 | passed |281.291 |
| Native update branch |14/14/0 | passed |206.872 |

These windows force an exact1 V endpoint; it must be shown separately from the prior full sweep's within-range point counts. Independent requalification confirms that AND and native_update have identical numerical IV rows and bit-identical final HDF datasets. Both have report-only carrier/continuity unsatisfied2/8; OR has8/10, with7 residual-only/2 update-only/7 both and one zero-update point. All original enforced gates pass, no ended-case solver remains, and each attempt+record+other partition closes. No uniform Newton improvement is established by switching stopping alone on this path. The two different original OR failures still need their trace/quality classification in the terminal analysis. Shared-VM times are not dedicated-machine speedups.

Cold0-100 V AND/OR controls are actually running; the native_update arm is queued under the two-worker limit. Each case has an independent14400-second safety cap. A controlled warm-target follow-up is also queued after these scans: use a qualified accepted state below each failed0-1 V OR target, then compare all three rules from exactly that same state/target. This is explicitly a controlled warm experiment, not a claimed replay of the original historical predictor.

## Required terminal analysis

Read first-pass iteration numbers from `stop_audit.jsonl`, keeping original stopping, native-update passing and independent physical acceptance separate. Exclude rejected/model-only candidates and terminal alpha0 from actual-writeback statistics. Compare packed states/raw directions and final states on fixed targets; compare adaptive regions[0,1],(1,5],(5,10],(10,100] with accepted points by actual outer voltage and retries/failures by target outer voltage. Do not infer direct failures by subtraction when a field is absent. Report exact endpoint clipping, OR branches and report-only conservation separately.

For every finished case, qualify solver/timeout/endpoint/original independent gates and absence of its own solver. Complete the cost partition identity attempt+record+other=total; inclusive substages are not additive and missing per-attempt times must stay missing. Necessary terminal evidence is collected once and member hashes verified. No success log is read; a failed preparation log is bounded-read only after exact terminal confirmation, with a read-once marker.

## Verified stage evidence and remaining work

The completed AV export,12 MP rows, frozen counterfactual,15 native-metric reconstructions and four R2 contract checks have exact successful supervisor/worker qualifications and no remaining stage processes. All6155 source files and the production binary/core/material identities were checked unchanged. The97-data-file stage archive is4120223 bytes, SHA256 `a8c055879a696ec6c62a700be97da36f724abb725cff8c6ca354c45dcfd1d08d`; archive and all members were verified. Successful logs were not read.

An additional AV-only frozen current-construction control has been launched: compute AV current directly from populations/potential with the already-existing population SG primitive, while leaving the original production transport, source-value qualification and solver dispatch intact. It is a constitutive arithmetic diagnostic, not a new production strategy. Its new binary/output directory does not overwrite the completed AV results. This control completed with finite fields and unchanged native nodes/production J: first-frame matrix difference remains1.53980233e-3 and next-iteration difference1.34788055e-7. Thus direct population SG arithmetic did not remove the six-row discrepancy either. Its26-file supplemental archive is1589474 bytes, SHA256 `fc939e83f5878cda4b1751c019a062ee07cf172870dd06627e9ed2cab752e6fe`, every member verified; exact export/analysis supervisors exited0 and the core stayed unchanged. The bounded AV experiment set is complete, with the native internal mechanism still not uniquely identified. The0-100 V scans, controlled warm targets, final shadow/quality/cost comparison and their terminal evidence remain pending. The native RHS branch remains unqualified and the native near-zero AV mechanism is not uniquely established. Monitoring continues every15 minutes in Chinese until the authorized work and honest qualification/limitations are recorded.


## Additional RHS physical-field audit

Existing NewtonPlot files contain PoissonRhs (C), eContinuityRhs/hContinuityRhs (A) and generic QuasiFermiRhs (A). A VM audit is now comparing these fields with the native algebraic matrix RHS. Initial filename-based matching used the wrong collection offset and was rejected by an explicit state check (up to0.0207 V disagreement); no normalization conclusion uses that rejected input. Corrected candidates cover12 frames, with3 source plots absent locally. Their physical potential differences from the frozen native matrix states range from about1e-15 to1.62e-10 V, so they are not claimed bit-identical. The audit retains units, phase uncertainty and per-block ratios; it does not grant an RHS stopping qualification. No new native solve is run to fill missing data.


## RHS authority/format audit update

The frozen parameter-audit receipt identifies the exact native pre/post filenames and recorded state-match bounds. Reading those authoritative files showed Ogawa format (rather than HDF5); the new field-reader preparation therefore exited1 before performing a scientific RHS comparison. Its exact supervisor exit was confirmed and its bounded failure read once. This is a format/preparation failure, not a numerical solver failure. The nearby local HDF candidates remain unqualified for exact near-root RHS matching.

A read-only probe of the installed U-2022.12 TDX converter was made in the isolated Sentaurus environment. Its initialization exited255 and reported no free TDX license. This does not establish that all Sentaurus products or licenses are unavailable; only this converter invocation failed to obtain access. No native solve or Vela simulation was started/restarted, and no running experiment was changed. Correctly decoded authoritative RHS fields remain pending; neither near-state medians nor fitted thresholds are promoted to the native RHSMin acceptance branch. The independently running update-branch controls continue.


## Qualified0-100 V controls and authoritative RHS decoding update

The original AND and OR bounded0-100 V controls completed without timeout, each with solver exit0, all accepted states passing its applicable original criterion and independent physical gates, and no remaining solver for that case. The shared observer worker is still running the native_update arm; no final batch worker exit is asserted. Successful logs were not read.

| Arm | Accept/attempt/direct fail | Solve time s | Final inner V | Report carrier/continuity unsatisfied |
|---|---:|---:|---:|---:|
| Original AND |123/165/42 |4407.730 |6.648838267668 |25/8 |
| Original OR |137/184/47 |4442.181 |6.648838267650 |55/10 |

Both end at external100 V and current approximately9.335116e-6 A/um. OR has89 residual-only,2 update-only and46 both branches, including one zero-update point. This is a bounded interval result, not complete-BV qualification.

| External interval V | AND accept/attempt/fail | OR accept/attempt/fail |
|---|---:|---:|
|[0,1] |13/13/0 |15/17/2 |
|(1,5] |12/16/4 |23/31/8 |
|(5,10] |25/35/10 |28/38/10 |
|(10,100] |73/101/28 |71/98/27 |

Each arm includes one exact100 V endpoint and one successful endpoint attempt. Top-level attempt+record+other closes exactly for both: AND4194196548891+207909261208+5500962004=4407606772103 ns; OR4205387394980+231305228642+5401469969=4442094093591 ns. Inclusive substages are not additive. Shared VM wall times are not exclusive speed comparisons. The32-data-file control archive is1877797 bytes, SHA256 `2eb29fcb5fea7c1e7447e4ac047af36467f1019a989f85502ef8606f77dbe76c`; archive and every member were verified.

The converter access issue above did not become a permanent evidence blocker. An independent read-only bounded Ogawa container reader and standard Zstandard scalar-array decoder recovered the authoritative existing TDR fields. Its container interpretation follows the public primary [Ogawa group implementation](https://github.com/alembic/alembic/blob/master/lib/Alembic/Ogawa/IGroup.cpp) and [data implementation](https://github.com/alembic/alembic/blob/master/lib/Alembic/Ogawa/IData.cpp); it does not invoke a vendor solver or alter an input. One isolated preparation revision referenced the wrong receipt name and was retained as failed before comparison; the corrected reader finished successfully.

All15 authoritative plots now pass explicit physical-state validation against the matching native algebraic exports. Maximum potential difference is1.33226763e-15 V and maximum carrier relative difference5.07145603e-15; these are serialization-level matches, not a claim of byte-identical independently serialized state. Region topology from the same-mesh prior HDF is checked by counts and full decoded state agreement; no near-state HDF values are substituted. Input plot hashes are retained.

Native physical PoissonRhs fields have unit C and the algebraic vector uses factor1/q=6.241449214576031e18, exactly matching the documented legacy q in these frames. Native e/hContinuityRhs fields have unit A and an observed native-only factor1.24325351991765e11; within-frame row-ratio spread is at most4.88498131e-15. The generic QuasiFermiRhs field is supported on11 drain aliases, and its summed contact relation is retained with cancellation sensitivity. These descriptive native field relations are not a calibration of Vela RHS, nor an acceptance threshold fitted to results. Agreement of Vela near-root physical residuals and all stopping branches remains unqualified, so RHSMin stays shadow-only. The26-data-file authoritative-field packet is6545839 bytes, SHA256 `b2fba451c3d3e3bbc649bed5560e11a6fcc97bcdf34f051279d44ea9cb5e673f`; all members and archive SHA were verified.

The native_update0-100 arm, the already-queued controlled warm pairs and final stopping/raw-writeback/quality/cost analysis remain pending. Earlier running snapshots in this report are historical. Production constants/core/materials, physical gates, budgets and original A1 admission block remain unchanged.


## RHS circuit-scale crosscheck and five accepted post states

A further read-only crosscheck separates a field's declared unit from its residual normalization. Simply treating a NewtonPlot continuity field labelled A as the unscaled one-micrometer terminal-current residual failed an explicit control and is not qualified. Under the frozen population/circuit coordinate convention, the raw native contact self-derivative10884.3831039366 divided by native VT gives voltage-row gain421026.8101838681. Combining that independent derivative with the specified1e7 ohm um resistor and production normalized-flux-to-current conversion2.98741982845799e-5 A/um gives carrier-row gain125778384.1055706. It agrees with the earlier independent-Jacobian calibration to2.22e-16 relative, without fitting RHS values.

The corresponding plotted-field/current-residual convention factor is33.86492002144134; the reexpression includes the documented RHS-versus-residual sign. This factor is neither a new physical constant nor a fitted stopping threshold, and contact-coordinate attribution remains specific to the frozen structure. For the12 previously exported hole rows, the corrected native algebraic absolute differences are at most3.354e-15 in the first frame and5.904e-15 in the next. Those local results do not qualify the complete RHS vector.

Existing compatibility-profile Windows residual exports at five accepted native post states were also reused, with no new solver invocation. Decoded authoritative post potentials and populations match those frozen initial-state values exactly in floating point. The comparison retains direct circuit-derived carrier/load scales and explicitly retains the prior independent-Jacobian Poisson scale.

| Accepted post state | Native plotted-field reconstructed RHS infinity norm | Vela mapped RHS infinity norm | Vector absolute difference infinity norm |
|---|---:|---:|---:|
|p01 |0.0249759 |0.0418370 |0.0402687 |
|p02 |0.00514560 |0.0420313 |0.0368857 |
|p03 |0.00225071 |0.0424745 |0.0402238 |
|p04 |0.00260627 |0.0425074 |0.0399012 |
|p05 |0.0248277 |0.0408630 |0.0349686 |

The largest vector difference is the circuit row in each sample. These native post norms are reconstructed from plotted fields using the pre-frame-validated field convention; they are not newly printed log values or independently exported post linear-system vectors. The comparison also spans the prior Windows compatibility observation and the active Linux stopping study. The absolute discrepancies exceed1e-5, so a general near-root native RHSMin equivalence remains unqualified. The observed update-error branch continues as the actual third stopping arm; the RHS branch remains shadow-only. No acceptance threshold was changed.

The34-data-file supplemental packet includes the five authoritative post plots, frozen Vela post residuals and input provenance,12 original row exports, native matrix/vector inputs and reproducible analysis. It is2583502 bytes, SHA256 `c591ad63aa4ddd76938200849e4859659297136111ce0c3632c90990ea79834a`; archive and every member SHA were verified. No completed job was polled or success log read to make this comparison.


## Final stopping comparison and qualification

All six adaptive cases,18 fixed controls and six controlled warm observations have ended. The stopping worker, warm worker and scientific-analysis worker each have exact supervisor exit0 and worker exit0, with no remaining task solver/stage process. All6155 source files, the source manifest, production binary/core/material and diagnostic wrapper identities were verified. No production defaults, floor,40-update budget,1% protection, models or thresholds were changed.

|0-100 V arm | Accept/attempt/direct fail | Sum of reported attempt iterations | Actual writebacks | Actual bootstrap/Full writebacks |
|---|---:|---:|---:|---:|
| Original AND |123/165/42 |1663 |1448 |423/1025 |
| Original OR |137/184/47 |1703 |1469 |498/971 |
| Observed native update branch |123/165/42 |1663 |1448 |423/1025 |

Actual writebacks require positive applied alpha and exclude rejected, model-review and linear-review rows. Full includes actual full and trial rows: AND/native have460 full and565 trial writebacks. All165/184/165 attempts match their consecutive outer/step/retry trace groups one to one. Reported iteration sums, trace counts and actual writebacks describe different things; none is substituted for a measured factorization count. Rejected/review counters overlap and are not additive.

AND and native_update have identical165-row attempt tables and identical1915-row Newton diagnostic tables. Their1915-row updates tables differ only in950 textual convergence_reason entries; all other columns agree. Physical IV and per-point iteration columns agree exactly and final HDF datasets are bit-identical. Thus the equal work counts are supported by matching paths, not just coincident totals. Native_update's elapsed4184.644 s versus AND4407.730 s does not establish an algorithmic speedup on this shared VM. Its failure classes are the same as AND:29 linear-gate/solve,11 iteration-budget and2 protected bootstrap prediction failures; OR has37/9/1.

The native_update arm's nonzero accepted states have valid update RMS below1, device residual below1e-8 and the unchanged independent gates. Its report-only carrier/continuity counts are25/8, retained. A direct call to the AND helper rejects its OR-mode metadata, so that helper is not an applicable qualification for this arm; this metadata check must not be mistaken for a numerical quality failure. The applicable update and physical gates were checked independently.

Within the six early frozen-state triples, initial packed/RHS hashes and last-recorded first-iteration raw-direction hashes agree across all three rules; stopping iteration counts are3,2,1,3,3,4 for p01_i01,p01_i02,p01_i03,p02_i01,p03_i01,p05_i01 respectively. Raw0 captures are the last recorded solve at iteration0, not a complete inventory of every rejected/retried candidate. On the AND adaptive path, no recorded earlier update-plus-load shadow pass preceded the original joint criterion; its zero-anchor exception is treated separately. Unselected shadow states did not receive every independent quality evaluation, so no hypothetical acceptance is inferred from the criterion alone. The underlying OR shadow on other paths can pass residual before update, but path changes prevent interpreting all such differences as removable iterations of the AND path.

### Controlled warm outcomes

| Target external V | AND | OR | Native update branch |
|---|---|---|---|
|0.599872743591 | linear qualification failure after6 reported iterations | same | same |
|0.944707732386 |21 updates, accepted |13 updates, accepted by original OR |21 updates, accepted |

Both triples have identical initial packed/RHS and last-recorded iteration0 raw hashes. The second target's OR stop has F6.86360834e-9, RMS270.078 and load-4.72558e-8 V; it correctly passes its residual branch and independent KCL/load checks. AND/native stop at F1.50606487e-9, RMS0.769696 and load-2.44728e-8 V. The generic checkpoint report's strict_gates field hardcodes AND, so its false value for the valid OR result is explicitly not used as OR rejection. Report-only carrier/continuity values were not separately exported by this checkpoint entry and remain unavailable; KCL and enforced system gates were checked. Three negative observations are preserved as scientific failures, despite all observer processes exiting0. These are controlled accepted-state starts below the target, not reconstructed historical predictors or proof that an earlier scan failure was rescued by changing stopping rules.

The native RHS branch remains shadow-only: its near-root absolute vector equivalence was not established, and the verified update branch must not be called complete Sentaurus semantics. The bounded AV checks establish AD consistency at the declared arithmetic anchors but do not uniquely identify the vendor's near-zero-current internal convention. No global AV sign reversal, smoothing, material adjustment or additional complete BV run follows from these findings.

The final866-data-file packet is66410042 bytes, SHA256 `cecb8682ca0d220e668dcc3e4c7b19f70c465997772b3f17986fbe54ecfda94a`. It contains necessary fixed/warm/adaptive terminal data, captures, audit ledgers, analysis, wrapper sources/configuration and input identities; archive and every member SHA were verified after one download. The previously qualified AND/OR100, AV and RHS packets are explicitly referenced by SHA rather than redownloaded. Top-level attempt+record+other partitions close for all adaptive arms. Per-attempt times remain unavailable, and included substages are not summed. Successful logs and numerical-case failure logs were not read for this final qualification.

The authorized bounded investigation is complete with these limitations recorded. No automatic compute queue remains; periodic vm-av monitoring is stopped after recording the terminal result. The original A1 complete-BV admission block is unchanged.
