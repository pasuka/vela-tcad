# Genius NPN BJT transport/source spatial characterization

These quantities are diagnostic characterization, not asserted acceptance gates.
Current-density comparison uses Vela's Sentaurus-style nodal quasi-Fermi-gradient reconstruction in A/cm^2.
Recombination integrals use the common triangular mesh and a 1 um out-of-plane depth.

| Quantity | Selected nodes | P95 log-magnitude error [decade] | Normalized error |
|---|---:|---:|---:|
| Electron current density | 5611 | 0.270332 | vector RMSE 0.563251 |
| Hole current density | 2121 | 1.46231 | vector RMSE 0.860804 |
| SRH recombination | 1807 | 0.306324 | L1 0.231876 |
| Auger recombination | 1627 | 0.903977 | L1 0.439981 |

Comparison status: **characterization_only**.
