# Genius NPN BJT SDevice current-semantics validation

## Scope

VBE=0.70 V and VCE=3.00 V on the exact 5611-node / 10940-triangle mesh. The original 2121-node current mask, 180 nodes above the 0.5-decade diagnostic limit, and the accepted SDevice state are retained.

## Frozen-state A/B

| Variant | State max delta | Tail P95 current delta (dec) | Tail max (dec) | Result |
|---|---:|---:|---:|---|
| ElementEdgeCurrent | 0.000e+00 | 0.000000 | 0.000000 | No plotted-current effect |
| hMobilityAveraging=ElementEdge | 0.000e+00 | 0.016375 | 0.023776 | Too small to explain the tail |

`ElementEdgeCurrent` leaves both node and element plotted current exactly unchanged on the loaded state. The mobility edge-average option changes the tail by at most 0.024 decade, far below the observed >=0.5-decade discrepancy.

## Native SDevice element-to-node reconstruction

The raw SDevice element-vector encoding is normalized by the fixed factor `1e-6/sqrt(2)`. An independently fitted scale on the main-current mask differs from this factor by less than 8e-5 relative. Area-weighting incident triangle currents gives:

| Evaluation set | Magnitude P95 (dec) | Normalized vector RMSE | Cosine |
|---|---:|---:|---:|
| Main current mask | 0.023944 | 0.021153 | 0.999776 |
| Original 180 tail nodes | 0.045367 | 0.221734 | 0.976374 |

## Vela SG edge-flux cross-check

Vela's directed SG line flux is divided by the dual-face length, recovered to one vector per triangle, then area-projected to nodes. This uses the same Vela flux data and changes only the representation support.

| Comparison | Median magnitude ratio | P95 magnitude error (dec) | Cosine |
|---|---:|---:|---:|
| Vela-recovered vs SDevice, 598 tail-adjacent cells | 1.067532 | 0.424222 | 0.971475 |
| Vela area-projected vs SDevice, 180 tail nodes | 1.268730 | 0.121906 | 0.999825 |
| Vela area-projected vs current Vela node field, 180 tail nodes | 22.058030 | 1.980122 | 0.599755 |
| Vela area-projected vs current Vela node field, 1700 dominant nodes | 0.999946 | 0.045585 | 0.997613 |

The same Vela conservative fluxes reproduce the SDevice tail after cell recovery and area projection, but are a median 22.06 times the current Vela node-vector representation. This is strong evidence that the large weak-current tail is primarily a support/reconstruction semantic difference, not a 22x transport-current physics difference.

## Convergence sensitivity

| Variant | Tail P95 delta (dec) | Tail max delta (dec) |
|---|---:|---:|
| -Extrapolate | 2.237939e-09 | 2.446988e-08 |
| Digits=8 + ExtendedPrecision(128) | 1.929793e-09 | 2.340297e-08 |

## Assessment and limitation

- The no-extrapolation solve is numerically identical to the accepted state, so continuation initialization is excluded as the tail cause.
- The element-current and Vela SG-flux cross-check upgrades the earlier hypothesis to strong evidence for differing node-vector construction semantics.
- This does not identify Synopsys' proprietary vertex reconstruction formula exactly; the strict oracle remains terminal current, KCL, and conservative control-surface flux rather than low-energy nodal-vector tails.
- Keep the existing 0.5-decade node gate as a diagnostic record. It should not alone fail the device while the dominant-current region and conservative flux checks pass.
