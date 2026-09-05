# Genius NPN BJT spatial-state acceptance

The asserted gate uses the exact common mesh at VBE=0.70 V and VCE=3.00 V.
Potential is checked on all nodes. Carrier-density decade errors are gated only where the Sentaurus reference density is at least 1e10 cm^-3.
Full-domain density statistics remain available as characterization and cannot override the registered gate.

| Field | Masked nodes | RMSE | P95 | Maximum | Pass |
|---|---:|---:|---:|---:|---:|
| Electrostatic potential [V] | 5611 | 0.00109652 | 0.0022881 | 0.0150609 | True |
| Electron density [decade] | 5196 | 0.0029782 | 0.0065252 | 0.023241 | True |
| Hole density [decade] | 2207 | 0.00693824 | 0.00776585 | 0.100667 | True |

Overall spatial-state pass: **True**

## Hole-density localization

The largest full-domain error is 2.87083 decade; the registered reference-significant mask reduces it to 0.100667 decade.

The highest errors are concentrated in reference-low-density nodes; see `hole_density_top_errors.csv` and the density-bin table in the JSON report.
