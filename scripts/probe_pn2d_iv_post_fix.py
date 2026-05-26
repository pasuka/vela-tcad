"""Quick verification of the residual-consistent ContactCurrent fix on pn2d IV.

Reads pre-existing simulation_iv configs in build/pn2d_tdr_tie_probe/vela,
re-runs the default Cathode and the Anode-contact variants, and prints
per-bias Vela vs Sentaurus reference plus the two-contact continuity sum.
"""
from __future__ import annotations

import csv
import bisect
import subprocess
import sys
from pathlib import Path

BUILD = Path("build/pn2d_tdr_tie_probe")
RUNNER = Path("build/vela_example_runner.exe")
VELA_DIR = BUILD / "vela"


def load_reference():
    ref_b, ref_i = [], []
    with (BUILD / "reference_curves" / "pn2d_iv_reference.csv").open() as f:
        for r in csv.DictReader(f):
            try:
                ref_b.append(float(r["bias_V"]))
                ref_i.append(float(r["current_total"]))
            except ValueError:
                pass
    pairs = sorted(zip(ref_b, ref_i))
    return [b for b, _ in pairs], [i for _, i in pairs]


def interp(ref_b, ref_i, b):
    j = bisect.bisect_left(ref_b, b)
    if j <= 0:
        return ref_i[0]
    if j >= len(ref_b):
        return ref_i[-1]
    b0, b1 = ref_b[j - 1], ref_b[j]
    i0, i1 = ref_i[j - 1], ref_i[j]
    return i0 if b0 == b1 else i0 + (b - b0) / (b1 - b0) * (i1 - i0)


def run(cfg_name: str):
    print(f"[run] {cfg_name}")
    subprocess.run(
        [str(RUNNER.resolve()), "--config", cfg_name],
        cwd=VELA_DIR,
        check=True,
    )


def read_iv(csv_path: Path) -> dict[float, float]:
    out: dict[float, float] = {}
    with csv_path.open() as f:
        for r in csv.DictReader(f):
            try:
                b = float(r["bias_V"])
            except ValueError:
                continue
            if r.get("handoff_stage", "") != "newton":
                continue
            out[round(b, 4)] = float(r["current_total_A_per_um"])
    return out


def main() -> int:
    if not RUNNER.exists():
        print(f"runner not found: {RUNNER}", file=sys.stderr)
        return 1
    ref_b, ref_i = load_reference()

    run("simulation_iv_default.json")
    run("simulation_iv_anode_contact.json")

    cath = read_iv(VELA_DIR / "pn2d_iv_default.csv")
    anode = read_iv(VELA_DIR / "pn2d_iv_anode_contact.csv")

    print()
    print(f"{'bias':>6} {'I_cath':>12} {'I_anode':>12} {'sum':>12} {'|I_cath/ref|':>12}")
    for b in sorted(cath):
        if b <= 1e-9:
            continue
        ic = cath[b]
        ia = anode.get(b, float("nan"))
        rf = interp(ref_b, ref_i, b)
        ra = abs(ic / rf) if rf else float("nan")
        print(f"{b:6.3f} {ic: .3e} {ia: .3e} {ic + ia: .3e} {ra:12.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
