# Genius NPN BJT hole-density discrepancy diagnosis

Bias: VBE=0.70 V, VCE=3.00 V; exact common mesh, no interpolation.

## QF decomposition

- Maximum hole-density error: 1.68359997 decade.
- The -99 to -101 mV driving-potential plateau contains 200 nodes; 200 also exceed 1 decade error.
- For SDevice p<1e2 cm^-3, the density-relation residual RMSE is 0.000294350074 decade.
- At the maximum-error node, measured and QF-predicted errors are -1.68359997 and -1.68359342 decade.

The low-density discrepancy is therefore explained by the difference in hQF-psi to numerical precision. This establishes the immediate state-variable cause; solver-limit causality requires the separate A/B runs.
