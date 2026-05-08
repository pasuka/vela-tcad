"""Run the PN diode DC sweep through the Vela Python API."""

from pathlib import Path

import vela


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    output_dir = repo_root / "examples" / "pn_diode" / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    sim = vela.run_dc_sweep(str(repo_root / "examples" / "pn_diode" / "simulation.json"))
    print(f"computed {len(sim)} DC sweep points")
