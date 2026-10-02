# SimpleMOS M9 chart map

| Report segment | Analytical question | Family / type | Fields | Supported takeaway | Delivery |
|---|---|---|---|---|---|
| Reference-density sensitivity | At what reference density does the Sentaurus response become material, and does comparison error improve? | Trend / multi-series line on logarithmic density | x: RefDens; y: maximum error change or maximum Sentaurus response; series: device and drain bias | `1e2 cm^-3` is effectively null; `1e6 cm^-3` and above worsens all four conditions | `simplemos_m9_refdens_error_change.png`, `simplemos_m9_refdens_sentaurus_response.png`, native report chart |
| Bias localization | At which gate biases does the strongest density control change current? | Trend / four-panel line | x: gate voltage; y: Sentaurus response and fixed-Vela residual; facet: device and drain bias | The strongest response is concentrated in deep-off and weak-inversion points | `simplemos_m9_refdens_pointwise_response.png` |
| Vela limiter strength | How strong is the production HFS limiter across the selected states? | Comparison / paired multi-series line | x: gate voltage; y: minimum edge mobility limiter; panel: carrier; series: device and drain bias | High drain bias creates much stronger local limiting than low drain bias | `simplemos_m9_vela_edge_limiter.png`, native report chart |

Palette policy: hard two-root cap within each individual comparison, with blue/gold/orange/pink used only to preserve the four fixed device/bias identities across the multi-series diagnostics. Marker shape and line style provide non-color distinction.
