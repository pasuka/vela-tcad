# Coarse-grid M1 comparison figures

The spatial maps use the common 5611-node/10940-triangle mesh at VBE=0.70 V
and VCE=3.00 V. Each field uses a shared SDevice/Vela color scale; the third
column is Vela minus SDevice. Current-density maps use the opt-in cell-first SG
node recovery and the registered SDevice magnitude mask of `1e-6` of peak.

The centerline figure samples the common triangular fields at `x=3.00 um`.
Terminal curves use all 31 collector-bias points from 0 to 3 V. The current
recovery P95 curve uses the four formal A/B biases at 0, 1, 2, and 3 V.

Reproduce from the worktree root:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\plot_genius_bjt_coarse_comparison_suite.py
```

`figure_manifest.json` records source/output SHA-256 hashes and spatial/vector
error metrics.
