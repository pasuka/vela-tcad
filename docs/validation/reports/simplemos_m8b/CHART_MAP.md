# SimpleMOS M8-B chart map

| Report segment | Question | Chart family / type | Encodings | Claim supported | Artifact and provenance |
|---|---|---|---|---|---|
| Sentaurus control response | Which HFS controls visibly change Sentaurus Id-Vg, and at which bias? | Small-multiple line chart | x: gate voltage; y: log10 current ratio to default; color: control; facet: device and drain voltage | Eparallel and RefDens produce material responses; contact and boundary controls are much smaller | `simplemos_m8b_sentaurus_response.png`; 20 direct Sentaurus curves, 51 exact Vg points each |
| Fixed-Vela error effect | Does a control consistently improve the Vela/Sentaurus comparison? | Grouped horizontal bar | y: device and drain voltage; x: change in maximum absolute log error; color: control | No tested control improves all four conditions; Eparallel and RefDens=1e8 worsen all four | `simplemos_m8b_error_change.png`; frozen M8-B comparison summary |
| Pointwise residual | Where in the gate sweep does each control change the residual? | Small-multiple line chart | x: gate voltage; y: absolute log current residual; color: control; facet: device and drain voltage | Large differences concentrate in low-current/transition regions; direction remains condition-dependent | `simplemos_m8b_residual.png`; fixed M8-A Vela curves versus M8-B Sentaurus curves |
| Best mean control Id-Vg | How large is the apparent best-mean improvement on the actual current curves? | Log-scale small-multiple line chart | x: gate voltage; y: absolute drain current; color/style: default Sentaurus, QF-at-contacts, fixed Vela | QF-at-contacts is only best by mean error change and is nearly coincident with default; it is not a uniform correction | `simplemos_m8b_best_control_idvg.png`; exact-bias current curves, no interpolation |

Color identity is held constant across the four figures: default/null gray, Eparallel blue, QF-at-contacts purple, no-parallel-boundary orange, and RefDens green.
