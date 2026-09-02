# Genius NPN BJT spatial-state acceptance

The asserted gate uses the exact common mesh at VBE=0.70 V and VCE=3.00 V.
Potential is checked on all nodes. Carrier-density decade errors are gated only where the Sentaurus reference density is at least 1e10 cm^-3.
Full-domain density statistics remain available as characterization and cannot override the registered gate.

| Field | Masked nodes | RMSE | P95 | Maximum | Pass |
|---|---:|---:|---:|---:|---:|
| Electrostatic potential [V] | 5611 | 0.00114132 | 0.00228107 | 0.0150803 | True |
| Electron density [decade] | 5196 | 0.00350584 | 0.00746252 | 0.0361364 | True |
| Hole density [decade] | 2207 | 0.0113454 | 0.0256416 | 0.100793 | True |

Overall spatial-state pass: **True**

## Hole-density localization

The largest full-domain error is 1.6836 decade; the registered reference-significant mask reduces it to 0.100793 decade.

The highest errors are concentrated in reference-low-density nodes; see `hole_density_top_errors.csv` and the density-bin table in the JSON report.
