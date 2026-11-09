# T470p Sentaurus U-2022.12 production A0: cold 0 to100 V

User requested AND and OR on T470p with emphasis on accepted points and target attempts by external-voltage interval. ProductionA0 automatic+protected model prediction review/rho reporting is retained; this is not always-Full or A1. The independent Vela VM full BV task is unchanged.

## Execution contract

Reuse the already qualified T470p compatibility binary/core and regenerated compatible materials, checking all source/input/binary hashes. No repeat build or regression/cold solve is scheduled: the identical previously qualified AND/OR cold0→.01V inputs and binary must be independently rechecked as frozen evidence before reuse. Material, binary and top-level physical_constants_profile must all be sentaurus_u2022_12. Two isolated cold0→100V windows run concurrently, each one PARDISO thread, with43200s safety limit and endpoint100V. Original physics/gates/N40/1%/floor1e-100 and Quasistationary controls remain; paths, profiling output and endpoint change only. This is not complete BV qualification and does not establish breakdown voltage.

## Current status

Connection recovered after user reboot and original strict host-key verification. Independent source/binary/core/input identities and frozen cold qualifications passed. Both actual windows started Beijing2026-11-03 14:00:52, deadline2026-11-04 02:00:52 each. Both windows have ended and terminal evidence is verified: AND failed before100 V; OR passed the 0–100 V window. The pair did not both succeed. See the final qualification below.

## Accounting

Use[0,1],(1,5],(5,10],(10,100]V without double-counting boundaries. Accepted points belong to their actual accepted external voltage; each completed target attempt, including retries/failures, belongs to its target external voltage. The zero-bias point/attempt is included in the first interval. Warm startup is not used. Counts exclude old cold preflight. Running counts describe completed rows only; do not infer failed attempts by subtraction. A region is fully covered only once its upper external voltage is reached. On early failure preserve last accepted and failed target and compare only common covered intervals; unvisited rows do not imply completed zero-cost regions.

The regional counter reproduces the previous repaired VM complete-BV terminal CSV exactly for bothAND/OR and passed explicit0/1/5/10/100 boundary tests. The comparison baseline is:

| Outer V | Native | Old AND | SRH-fixed AND | Old auto OR | SRH-fixed auto OR |
|---|---:|---:|---:|---:|---:|
|0–1|13/13|18/22|16/18|19/23|21/27|
|1–5|11/13|21/27|15/19|20/26|14/16|
|5–10|8/9|14/17|43/59|29/40|22/30|
|10–100|12/12|71/98|74/103|71/99|71/99|

Cells are accepted/attempts; failures are reported separately for the new task. Terminal reporting will add compatibleT470pAND/OR columns and counts with/without the forced exact100V endpoint. The older complete sweeps may cross100V without accepting exactly100V; endpoint clipping must be separated from policy workload differences.

## Interpretation and qualification

A0AND uses enforce_joint, OR uses enforce, both retain load enforce/KCL/nonfinite/all original independent gates; report-only carrier/continuity failures and OR residual-only/update-only/both/zero-update branches remain explicit. Success of this window requires solver0/no timeout/endpoint100V/all accepted applicable gates and no owned solver residual process. A failed window retains evidence and counts and will not be silently restarted or given relaxed tolerances/budget. Success logs will not be read; failure logs only if needed, bounded once after termination. Terminal scientific evidence is archived once and each file SHA checked.

The comparison combines a new T470p run with older VM/native evidence. Compiler/platform/roundoff and endpoint differences are retained as limitations; counts alone cannot identify a unique constants/SRH cause or dedicated-machine speedup. No new Sentaurus simulation or extra CODATA control run is authorized here.

## Connection restored after user reboot

The saved computer-name.local route recovered and strict existing SSH host-key verification passed. R1 WMI startup returned invalid-parameter before any worker was created; that receipt is preserved. R2 changed only the independent launcher to native Windows paths and hidden startup, then successfully created the new worker. The first startup snapshot was verify_reused_identity; subsequent actual window start/deadline is recorded below. The task is not being rebuilt or restarted, and the VM full scans continue independently. Detailed receipts remain ignored.

## Actual windows confirmed

Both identity/preflight-reuse checks passed. AND andOR actually started Beijing2026-11-03 14:00:52, each deadline2026-11-04 02:00:52 (43200s); prior source/input verification is outside this solver budget. At about80s AND had9accepted/10attempts/0directfailed, lastouter0.3132512447185/inner0.2893941688225V/I2.3857079358e-9A/um; OR10/13/1, last0.437833177002/0.413023894612V/I2.480926012e-9A/um. The0–1V region remains incomplete; higher ranges are unvisited. IV may lag diagnostic attempt records, so differences between these counters must not be interpreted as failures. No window success or count-reduction claim yet.

At about220s, AND10accepted/14attempts/3failed through outer0.4253962176065V (0-1incomplete). OR19/27/6 through1.244190022121V; its completed0-1region is15/17/2, compared with the earlier repaired VM21/27/6. OR1-5is partial4/10/4. These are partial-window records with cross-platform limitations, not a complete0-100V result or unique constants-causation claim.


## Monitoring snapshot, 2026-11-03 about14:10 Beijing

Both0-100V windows remain running, elapsed about541s, deadlineNov4 02:00:52. AND32accepted/40attempts/7failed, lastouter3.148644071160/inner3.077848019473V/I7.079612538e-9A/um. OR37/49/11, last2.797474768796/2.735717309440V/I6.175756463e-9A/um. Completed0-1V region counts (accepted/attempts/failures):

| Mode | Compatible T470p | Previous SRH-fixed CODATA VM |
|---|---:|---:|
|AND|24/29/5|16/18/2|
|OR|15/17/2|21/27/6|

AND has more and OR fewer points/attempts in this completed region; there is no uniform reduction. The1-5V region is partial: AND8/11/2, OR22/32/9. Higher regions are unvisited. Cross-platform/compiler/roundoff limitations remain; no0-100V success or unique constants causation is claimed.


## AND terminal and requested failure explanation

AND finished numerically unsuccessful, not timed out:41accepted/66attempts/25failed, solver1/group worker4. Last accepted outer5.580805160826V/inner5.312860934612V/I2.679442456146e-8A/um. The final target5.580820165011V, retry13, step15.0041848457uV, exhausted40iterations. F1.8750160313e-9 and |L|2.1256835403e-8V passed, but valid RMS156.273990921 failed the original AND requirement. Halving would produce7.50209uV belowMin10uV. All41accepted points independently pass originalAND/independent gates; report-only carrier/continuity unsatisfied counts16/6 retained. No ownedANDsolver remains; OR continues, so no batch final worker code is claimed. Terminal15-file archive6407a84fc8cbb01d30e803491e1bfb5aeddc84fb1b02998f3bc73224fb307b44 was downloaded once and every member SHA checked; the terminal failure log was bounded-read once only.

For the completed0-1V region, AND's five failures are all40iteration exhaustion. At targets0.5375411905,0.5179345111,0.7966027072,0.8097372117V, F andload passed but RMS was17.6923,491.0992,184.6004,26.5720 respectively. At0.5544003181V, F8.30261655e-8 andRMS50705.5447 both failed (load passed). No linear failure is inferred from terminal NaN fields.

OR's two0-1V failures were at0.5998727436V(iter6) and0.9447077324V(iter5). Bootstrap failed and same-stateFull retry also failed the scaled absolute linear-defect requirement:1.02056305e-4 and1.06723380e-3 versus1e-8. Componentwise errors8.70218e-17 and4.41531e-12 passed1e-10, but both linear gates are required. These were rejected nonlinear trial states, not accepted states.

The native13accepted/13attempts/0failure comparator is archivedT-2022.03-SP2, while the15same-state frames/constants audit useU-2022.12. Its12nonzero0-1Vtargets converged in3-8iterations; none of the seven rejected Vela targets is exactly one of those native targets. Native zero failures establish success on its own path and criteria; they do not establish acceptance of these same Vela trial states. Constants compatibility does not by itself equate effective Jacobian, protection/history, Full candidate handling, linear qualification or stopping semantics. The first divergence mechanism for these exact seven attempts remains unproven. No new Sentaurus run, relaxation or solver change was made for this analysis. Frozen low01 classification and original-target context are preserved in localJSON evidence.


## Regional progress, 2026-11-03 about14:31 Beijing

AND remains the independently qualified numerical failure described above and is not polled individually. OR remains running at1822.645s,120accepted/161attempts/40failed, lastouter54.39993009601/inner6.654572455732V/I4.774535748e-6A/um. Its final10-100Vregion is still incomplete. Counts are accepted/attempts/failures:

| Outer V | Old fixed VM AND | Compatibility T470p AND | Old fixed VM OR | Compatibility T470p OR |
|---|---:|---:|---:|---:|
|0-1|16/18/2|24/29/5|21/27/6|15/17/2|
|1-5|15/19/4|13/16/3|14/16/2|32/42/10|
|5-10|43/59/16|4/21/17(partial; stopped5.580805V)|22/30/8|23/31/8|
|10-100|74/103/29|not reached|71/99/28|50/71/20(partial)|

The completed0-5V totals are newAND37/45/8 versus old31/37/6, and newOR47/59/12 versus old35/43/8. Both totals increased despite the OR0-1Vreduction. These remain cross-host comparisons, not unique constants causation. Uncovered ranges are not complete zero-cost regions, and no0-100Vsuccess is claimed.

## Final qualification, 2026-11-03

OR ended at 2026-11-03T14:36:12.106129+08:00 Beijing after2119.815825 s (35min19.816s), not timed out, solver0/groupworker0. It reached exactly100 V with142accepted/191attempts/49failed; last internal voltage6.648838267773 V and current9.33511618534e-6 A/um. Every accepted point independently passed the original OR and independent gates. Branches:100 residual-only,2 update-only,40 both,1 zero-update point. Report-only carrier102/continuity9 unsatisfied points remain separate. At the final point F1.77794834e-9 andload1.21202305e-7 V passed; RMS4.65485 alone does not pass AND and is not used to disqualify the original OR rule.

AND remains the previously verified numerical failure41/66/25 before100 V. The batch worker exited4 as the expected aggregate of failedAND and successfulOR; the worker and all task-owned Python/solver processes are absent. All6154 source files, binary/core and input hashes were checked again and unchanged. The OR terminal packet contains19 data files,657647 bytes, SHA256 `6aa595091b5f71896332f034bde940b2bad81ca32b233d0e9e91d6f33b314396`; archive and every member SHA passed. No successful solver log was read. Both T470p groups are now removed from polling; the separate VM full sweeps continue.

Final regional counts are accepted/attempts/direct failures:

| Outer V | Old SRH-fixed VM AND | Compatible T470p AND | Old SRH-fixed VM OR | Compatible T470p OR |
|---|---:|---:|---:|---:|
|0–1|16/18/2|24/29/5|21/27/6|15/17/2|
|1–5|15/19/4|13/16/3|14/16/2|32/42/10|
|5–10|43/59/16|4/21/17; partial, stopped5.580805 V|22/30/8|23/31/8|
|10–100|74/103/29|not reached|71/99/28|72/101/29|

OR totals over0–100 V are142/191/49 versus old128/172/44:14 more accepted points,19 more attempts and5 more failures. Its forced exactly100 V endpoint contributes1 accepted point/1 successful attempt; after removing that endpoint the counts are141/190/49, still13/18/5 more than the old regional totals. This exclusion is accounting only, not a replay of the un-clipped path. The endpoint-clipped10–100 V row becomes71/100/29. AND only has fully comparable completed0–5 V intervals (37/45/8 versus old31/37/6); its incomplete5–10 V region cannot support full-region workload conclusions.

OR top-level inclusive cost closes: attempt1975.7876897 s + recording138.5270542 s + other4.3149559 s = total2118.6296998 s. Solver-process wall time also includes launcher/monitor overhead. Nested counters must not be added to these top-level values, and per-attempt elapsed times are unavailable. Raw qualification, config, IV/updates/attempts/Newton, profile/material and performance evidence are archived once; the previously missing failed-AND partition remains missing.

There is no uniform reduction in work: OR0–1 V improves its counts while1–5 V worsens; totalOR0–100 V costs more attempts, andAND did not finish. Cross-host/compiler/roundoff and exact-endpoint differences preclude unique constants/SRH attribution. This is a 100 V external-voltage window result, not complete breakdown qualification or a new Sentaurus sweep. The native regional reference remains archived T-2022.03-SP2, distinct from the U-2022.12 constants/frozen-frame audit.
