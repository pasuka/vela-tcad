# SimpleMOS M46 full-matrix requalification

The complete eight-device, sixteen-curve, 816-point production SimpleMOS matrix was rerun after M43/M44. All curves passed and every portable comparison CSV is bitwise identical to the frozen M8 result.

- Maximum cross-TCAD error: 0.109418680924 dex
- Maximum relative error: 0.28652633602
- Maximum endpoint error: 0.0131744123411 dex
- Maximum current shift versus M8: 0 dex
- Worst point: n23_vd_0p05, Vg=0.05 V

M43/M44 therefore close the identified diagnostic and frozen-replay inconsistencies without changing the default BGN-on Id-Vg matrix. The remaining 0.1094 dex cross-TCAD peak is unchanged and requires a different causal target.
