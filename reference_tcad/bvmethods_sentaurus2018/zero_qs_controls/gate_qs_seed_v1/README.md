# Auxiliary gate QS acquisition candidate

This separately authorized native Sentaurus control starts at zero bias,
ramps the drain outer voltage to 13000 V, sweeps the gate 0 -> 0.2 -> 0 V,
then resumes the original drain voltage goal of 100000 V. All stages use
Quasistationary with coupled Poisson/electron/hole equations. There is no
transient acquisition or loaded seed.

The model declarations, material parameters, fixed `1e7 ohm*um` resistor and
`1.443e-3 A/um` stop remain the current repository baseline. To reconstruct
the input directory, use `source/models.par` from the BVmethods parent and
the byte-verified original mesh recorded in the linked report. Generated
TDR/PLT/log/Save files remain on the VM and in ignored build mirrors.

This changes the initialization bias path. Its qualification is separate
from the original fixed-gate, zero-origin voltage QS trajectory. Native
acquisition does not establish dynamic stability or Vela/native field
agreement. The source SHA-256 is
`5a7d5468a468ff82538906fb147dff854373ae4e49812fb924ee2188dd608283`.

Independent checker: `scripts/validate_wp14_gate_qs.py`, registered as
`wp14_gate_qs`. Results and archive receipt:
[zero-QS diagnosis](../../../../docs/validation/bvmethods_zero_qs_diagnosis_2026-10-06.md).

Executed on 2026-10-06: failed. The gate ramp reached only 0.0966954759163124 V;
after restoring zero gate bias, the final drain ramp stopped at
0.00140926509666232 A/um below the declared threshold. All 203 accepted states
passed the independent physical closure checks; stage completion did not.
Raw data are sealed in `bvmethods_zero_qs_controls/v2`, not promoted into the
current v4 reference manifest.
