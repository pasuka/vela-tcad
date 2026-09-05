# Genius NPN BJT low-current hole-current tail diagnosis

## Scope

VBE=0.70 V and VCE=3.00 V, on the exact 5611-node common mesh. The original SDevice-relative 1e-6 peak mask and 0.5-decade node limit are retained.

## Main result

The accepted mask contains 2121 nodes; 180 exceed 0.5 decade. The failure is reported by strata below; an SDevice directed-edge flux is not available, so interpolation semantics cannot be eliminated as a strict oracle-level possibility.

## Reference-current magnitude

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| [1e-2,1] | 1700 | 0 | 0.000 | 0.0427 | 0.998 | 1.000 |
| [1e-3,1e-2) | 218 | 38 | 0.174 | 1.5010 | 0.002 | 0.000 |
| [1e-4,1e-3) | 88 | 39 | 0.443 | 1.6544 | 0.000 | 0.000 |
| [1e-5,1e-4) | 63 | 53 | 0.841 | 2.1297 | 0.000 | 0.000 |
| [1e-6,1e-5) | 52 | 50 | 0.962 | 1.8455 | 0.000 | 0.000 |

## Electrical region

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| n_emitter_side | 640 | 0 | 0.000 | 0.0174 | 0.018 | 0.024 |
| p_base | 1481 | 180 | 0.122 | 1.4740 | 0.982 | 0.976 |

## Depth band

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| [0,0.25) | 793 | 17 | 0.021 | 0.0603 | 0.681 | 0.246 |
| [0.25,0.50) | 886 | 19 | 0.021 | 0.0756 | 0.275 | 0.582 |
| [0.50,0.75) | 442 | 144 | 0.326 | 1.5923 | 0.043 | 0.172 |

## Lateral zone

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| base_contact_window | 338 | 21 | 0.062 | 0.6380 | 0.367 | 0.051 |
| base_emitter_gap | 429 | 32 | 0.075 | 0.5522 | 0.502 | 0.472 |
| left_of_base_contact | 111 | 12 | 0.108 | 0.8870 | 0.000 | 0.000 |
| right_base_extension | 318 | 21 | 0.066 | 0.6747 | 0.010 | 0.010 |
| right_of_base_profile | 122 | 30 | 0.246 | 1.7175 | 0.000 | 0.000 |
| under_emitter_window | 803 | 64 | 0.080 | 0.5898 | 0.121 | 0.466 |

## Boundary type

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| contact_base | 18 | 0 | 0.000 | 0.0888 | 0.362 | 0.020 |
| contact_emitter | 34 | 0 | 0.000 | 0.0502 | 0.007 | 0.001 |
| contact_one_ring | 52 | 0 | 0.000 | 0.0158 | 0.008 | 0.017 |
| insulating_exterior | 49 | 1 | 0.020 | 0.1050 | 0.276 | 0.047 |
| interior | 1968 | 179 | 0.091 | 1.0191 | 0.347 | 0.916 |

## Reference-current direction

| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |
|---|---:|---:|---:|---:|---:|---:|
| +x | 1004 | 141 | 0.140 | 1.4957 | 0.571 | 0.895 |
| +y | 280 | 12 | 0.043 | 0.4242 | 0.406 | 0.078 |
| -x | 115 | 0 | 0.000 | 0.0419 | 0.004 | 0.003 |
| -y | 722 | 27 | 0.037 | 0.1387 | 0.019 | 0.024 |

## Hypothesis checks

- SDevice nodal-vector edge-projection/recovery round trip: P95 0.0461 decade.
- Vela versus original SDevice: P95 0.8280 decade.
- Vela versus round-tripped SDevice: P95 1.0160 decade.
- Pearson correlation between signed current log-ratio and the common mu-p-grad(phi_p) proxy log-ratio: -0.15681585818466778.
- Pearson correlations of absolute current error with SDevice round-trip, hole density, mobility, and QF-gradient errors: 0.29592031145030784, 0.44835037226679036, -0.13233016599363004, and 0.08071988221149266.
- On the 180 tail nodes, the median SDevice current/local-proxy log ratio is 1.3472 decade; the Vela value is 0.0001 decade.

## Assessment

- Verified localization: all 180 tail nodes have lower Vela current; they are in the p-type base, and 179 are interior nodes. The 0.50-0.75 um depth band contains 144 tail nodes.
- Recovery-weight hypothesis is not supported: uniform, primal-length, dual-face, dual-face-squared, and dual-area fits span only 0.0069 decade in P95.
- Local-state/model-factor mismatch is not supported as the primary explanation: on tail nodes, the P95 absolute differences are 0.0322 decade for p, 0.0043 decade for mobility, and 0.0070 decade for the common-mesh QF gradient.
- The strongest evidence points to SDevice's weak-current nodal field construction: its tail-node current is a median 22.2 times the local q*mu*p*|grad(phi_p)| proxy, while Vela is 1.000 times. This is a likely interpretation, not a proof, because SDevice's directed edge flux is unavailable.

## Data quality and limitation

All SDevice fields contain one finite, ordered row for every common-mesh node. The Vela state/VTK and SDevice field-manifest hashes are recorded in the JSON result. The SDevice current is an exported nodal vector, while Vela's authority is a conservative SG line flux; therefore this analysis distinguishes evidence for the two hypotheses but does not claim access to an unavailable SDevice edge-flux oracle.
