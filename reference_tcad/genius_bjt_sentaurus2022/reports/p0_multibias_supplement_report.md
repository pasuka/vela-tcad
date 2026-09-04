# Genius NPN BJT P0 multi-bias supplement

## Technical summary

The missing VCE=2 V SDevice state was regenerated with Sentaurus Device
T-2022.03-SP2 from the accepted 5611-node mesh and the same M1 physics deck.
The four audited states now cover VCE=0, 1, 2, and 3 V at VBE=0.70 V.

The added points support the existing diagnosis.  Vela's conservative SG
transport remains consistent with the SDevice terminal solution: the absolute
Ic, Ib, and Ie errors are below 0.8% at all four biases, while the Vela terminal
KCL defect is between 1.3e-20 and 2.1e-16 A/um.  The registered spatial-state
masks also remain inside their 3 V thresholds at every audited bias.

The only metric that degrades strongly with VCE is the logarithmic magnitude
tail of the exported nodal hole-current vector.  Its P95 error increases from
0.016 decade at 0 V to 0.828 decade at 3 V, even though vector NRMSE stays near
0.062 and cosine similarity stays above 0.9978.  Nodes above 0.5 decade carry
at most 3.93e-5 of the reference-current energy.  This separates a weak-current
node-representation discrepancy from a material conservative-flux error.

## Terminal currents remain within 0.8%

| VCE (V) | Ic abs. error | Ib abs. error | Ie abs. error | Vela KCL abs. (A/um) |
|---:|---:|---:|---:|---:|
| 0 | 0.291% | 0.322% | 0.475% | 1.39e-20 |
| 1 | 0.788% | 0.480% | 0.781% | 2.02e-16 |
| 2 | 0.787% | 0.482% | 0.780% | 1.27e-17 |
| 3 | 0.785% | 0.483% | 0.778% | 1.00e-16 |

These terminal currents are the observable conservative-flux comparison.  A
native SDevice directed-edge flux is not exported, so an edge-by-edge SG oracle
is unavailable.

## Spatial state remains close on the registered masks

Density errors below use log10 decades and the registered SDevice-reference
mask of at least 1e10 cm^-3.  Potential uses all common nodes.

| VCE (V) | Potential RMSE (V) | Potential max (V) | Electron P95 (dec) | Electron max (dec) | Hole P95 (dec) | Hole max (dec) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 4.90e-5 | 1.08e-3 | 0.00117 | 0.0169 | 0.00222 | 0.0190 |
| 1 | 6.00e-4 | 7.81e-3 | 0.00796 | 0.0414 | 0.0169 | 0.1312 |
| 2 | 8.75e-4 | 1.20e-2 | 0.00681 | 0.0347 | 0.0107 | 0.2019 |
| 3 | 1.10e-3 | 1.51e-2 | 0.00653 | 0.0232 | 0.00777 | 0.1007 |

The 0, 1, and 2 V rows are characterization against the existing 3 V contract;
they do not register new production gates.  All would pass the same numeric
limits.

## Hole-current tail grows, but its energy contribution stays negligible

The current-vector mask retains SDevice nodes whose magnitude is at least
1e-6 of the per-bias peak.

| VCE (V) | Selected nodes | Hole-current P95 (dec) | Vector NRMSE | Cosine similarity | Nodes >0.5 dec | Reference-energy share of tail |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5611 | 0.0164 | 0.0658 | 0.997893 | 0 | 0 |
| 1 | 2492 | 0.2474 | 0.0624 | 0.998085 | 23 | 2.29e-9 |
| 2 | 2257 | 0.6167 | 0.0622 | 0.998103 | 121 | 3.92e-5 |
| 3 | 2121 | 0.8280 | 0.0622 | 0.998104 | 180 | 3.63e-6 |

The C++ dual-control-face reconstruction agrees with an independent Python
replay to a maximum of 4.84e-13 A/cm2 across the four biases.  Thus the observed
tail is not an implementation mismatch between the documented reconstruction
and its C++ output.  In the base-collector window, the recovered nodal vector's
NRMSE against the SDevice nodal vector grows moderately from 0.0256 to 0.0628.
By contrast, directly comparing a raw SG edge projection with a SDevice nodal
projection grows from 0.0111 to 0.7189.  The latter compares different supports
and must not be interpreted as a conservative-flux validation metric.

## Recombination differences are stable above 1 V

| VCE (V) | SRH integral Vela/SDevice | SRH normalized L1 | Auger integral Vela/SDevice | Auger normalized L1 |
|---:|---:|---:|---:|---:|
| 0 | 0.8596 | 0.1416 | 0.9700 | 0.0300 |
| 1 | 0.7788 | 0.2213 | 0.9878 | 0.0122 |
| 2 | 0.7787 | 0.2213 | 0.9878 | 0.0122 |
| 3 | 0.7787 | 0.2214 | 0.9878 | 0.0122 |

The SRH integral remains about 22.1% lower in Vela after 1 V and is a genuine
model-alignment item.  The Auger integral is already within about 1.2%.  Neither
trend tracks the rapidly increasing nodal hole-current log-tail, so neither is
a sufficient explanation for that tail.

## Method and data lineage

- SDevice states were verified from TDR-internal contact voltages, not inferred
  from filenames.  Every audited state has 5611 nodes and 10940 elements.
- Vela states are the P0 fixed-bias accepted states.  Current and recombination
  fields were regenerated only from those states.
- Spatial density statistics exclude SDevice reference density below
  1e10 cm^-3.  Current-vector statistics use a per-bias 1e-6-of-peak mask.
- The 2 V SDevice TDR SHA-256 is
  `8f1e876838531415311315825ad56c6637e5a40d41fce2aafac698566ed4134a`.

## Limitations and next step

This analysis cannot prove that the remaining tail is solely an SDevice nodal
interpolation effect because SDevice does not expose the internal directed-edge
flux used by its discretization.  It does show that the terminal conservative
observable is accurate, the accepted carrier state is accurate in significant
regions, and the tail has negligible reference-current energy.

The next engineering step should therefore keep the 0.5-decade tail metric as a
diagnostic, not an overall veto, while adding a registered conservative section
flux comparison.  Separately, the stable SRH integral offset should be pursued
through lifetime, doping dependence, and quadrature alignment.

