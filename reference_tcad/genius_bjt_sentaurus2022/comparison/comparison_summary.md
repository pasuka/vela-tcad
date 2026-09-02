# Genius NPN BJT Sentaurus/Vela comparison

All 31 requested collector biases use directly reported collector, base, and emitter currents.
Operational pass means exact bias selection, solver convergence, and terminal KCL passed.
M1 numerical parity is evaluated over VCE=0.5-3.0 V against the pre-registered maximum log-error thresholds.

| Model | Operational | Numerical parity | Sentaurus Ic @ 3 V (A/um) | Vela Ic @ 3 V (A/um) | Ic ratio | Ib ratio | Sentaurus beta | Vela beta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| M0 | True | None | 1.843315039e-06 | 3.397337261e-06 | 1.84306 | 0.444397 | 114.928 | 476.643 |
| M1 | True | True | 2.807532434e-06 | 2.790668781e-06 | 0.993993 | 1.05219 | 46.6464 | 44.0665 |

## M1 numerical parity gate

| Observable | Maximum allowed absolute log10 error | Observed maximum | Pass |
|---|---:|---:|---:|
| Ic | 0.05 | 0.00269368076 | True |
| Ib | 0.05 | 0.0220922063 | True |
| Ie | 0.05 | 0.0021204468 | True |
| beta | 0.05 | 0.0247346089 | True |

Overall operational pass: **True**
Asserted numerical parity pass: **True**
Asserted spatial-state pass: **True**
Overall pass: **True**
