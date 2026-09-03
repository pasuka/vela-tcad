# Genius NPN BJT transport/source spatial acceptance

The initial engineering gates below were registered after the original characterization run; they are regression gates, not an independent blind validation.
Current-density comparison uses Vela's Sentaurus-style nodal quasi-Fermi-gradient reconstruction in A/cm^2.
Recombination integrals use the common triangular mesh and a 1 um out-of-plane depth.

| Quantity | Selected nodes | P95 log-magnitude error [decade] | Normalized error |
|---|---:|---:|---:|
| Electron current density | 5611 | 0.270794 | vector RMSE 0.564365 |
| Hole current density | 2121 | 1.48371 | vector RMSE 0.815032 |
| SRH recombination | 1807 | 0.289153 | L1 0.2214 |
| Auger recombination | 1627 | 0.903665 | L1 0.424581 |

Electron current-density pass: **True**.
Hole current-density pass: **False**.
SRH recombination pass: **True**.
Auger recombination pass: **False**.
Comparison status: **asserted**.
Overall transport/source pass: **False**.
