# Genius NPN BJT Sentaurus/Vela comparison

All 31 requested collector biases use directly reported collector, base, and emitter currents.
Operational pass means exact bias selection, solver convergence, and terminal KCL passed.
Numerical parity is intentionally characterization-only because WP0-WP2 did not pre-register a Vela error tolerance.

| Model | Operational | Sentaurus Ic @ 3 V (A/um) | Vela Ic @ 3 V (A/um) | Ic ratio | Sentaurus beta | Vela beta |
|---|---:|---:|---:|---:|---:|---:|
| M0 | True | 1.843315039e-06 | 3.397337261e-06 | 1.84306 | 114.928 | 476.643 |
| M1 | True | 2.807532434e-06 | 2.679005864e-06 | 0.954221 | 46.6464 | 9.77105 |

Overall operational pass: **True**
