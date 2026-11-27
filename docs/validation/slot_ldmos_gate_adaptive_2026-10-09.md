# SLOT adaptive gate ramp, 2026-10-09

Scope: new cold equilibrium, Vd=0, gate 0 to 5 V only. No drain scan.

Policy: first 0.005 V; Increment=1.5; reference iterations=25; Decrement=2; MaxStep=0.25 V; MinStep=5e-8 V; at most 100 attempts and 8 consecutive failures. Growth uses full Newton accepted update count (zero clamped to one); Poisson preparation iterations are separately recorded. This differs from native monolithic initialization. Newton limits and physical gates unchanged.

AverageBox, Fermi, full target physics and contact avalanche fallback; eps_row=.001, closure=.01, abstol=1e-9/reltol=1e-8. Frozen solver copied from s08averagecold.

Result: {"completed": true, "goal_reached": false, "stop_reason": "ValueError: Diagnostic predictor permits at most 0.1 V bias change", "elapsed_s": 930.1042113999999, "resources_unchanged": true}

Accepted counts: {"equilibrium": 1, "gate": 8, "drain": 0}

Record audit passed: True. This is not an independent numerical reproduction.

Evidence: build/s09gateadaptive/run/summary.json, preflight.json, audit.json, gate_steps.csv. Failed and preparation states are excluded.
