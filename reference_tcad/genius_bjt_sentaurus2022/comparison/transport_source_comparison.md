# Genius NPN BJT transport/source spatial acceptance

The initial engineering gates below were registered after the original characterization run; they are regression gates, not an independent blind validation.
Current-density comparison uses the dual-face-length-weighted nodal representation of Vela's production SG line fluxes in A/cm^2.
The SG line fluxes remain the conservative finite-volume authority; the node vectors are field-comparison diagnostics.
Recombination integrals use the common triangular mesh and a 1 um out-of-plane depth.

| Quantity | Selected nodes | P95 log-magnitude error [decade] | Normalized error |
|---|---:|---:|---:|
| Electron current density | 5611 | 0.0157728 | vector RMSE 0.0351488 |
| Hole current density | 2121 | 0.827987 | vector RMSE 0.0622053 |
| SRH recombination | 1807 | 0.289153 | L1 0.2214 |
| Auger recombination | 1627 | 0.903356 | L1 0.01219 |

Electron current-density pass: **True**.
Hole current-density pass: **False**.
SRH recombination pass: **True**.
Auger recombination pass: **True**.
Comparison status: **asserted**.
Overall transport/source pass: **False**.
