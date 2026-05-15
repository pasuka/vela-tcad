"""Run documented device curve examples through the Vela Python API.

Build with ``-DVELA_ENABLE_PYTHON=ON`` and set ``PYTHONPATH`` to the generated
package directory, for example ``PYTHONPATH=build/python/Debug``.
"""

from pathlib import Path

import vela


def _print_curve(label: str, points: list[dict[str, object]]) -> None:
    print(f"{label}: {len(points)} points")
    if not points:
        return
    first = points[0]
    print(
        "  "
        f"type={first['curve_type']} "
        f"bias={first['bias_contact']} "
        f"current={first['current_contact']} "
        f"csv={first['output_csv']}"
    )
    last = points[-1]
    if first["curve_type"] == "cv_quasistatic":
        print(f"  final capacitance={last['capacitance']} F/m")
    elif first["curve_type"] == "bv_reverse":
        print(
            "  "
            f"breakdown={last['breakdown_detected']} "
            f"voltage={last['breakdown_voltage']} V "
            f"criterion={last['breakdown_criterion']}"
        )
    else:
        print(f"  final total current={last['total_current']} A/m")


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    for device in ("pn_diode", "nmos2d_dd"):
        (repo_root / "examples" / device / "outputs").mkdir(parents=True, exist_ok=True)

    curves = [
        ("PN IV", vela.run_iv_curve(repo_root / "examples" / "pn_diode" / "simulation_iv.json")),
        ("PN CV", vela.run_cv_curve(repo_root / "examples" / "pn_diode" / "simulation_cv.json")),
        ("PN BV", vela.run_bv_curve(repo_root / "examples" / "pn_diode" / "simulation_bv.json")),
        (
            "NMOS Id-Vd",
            vela.run_iv_curve(repo_root / "examples" / "nmos2d_dd" / "simulation_iv.json"),
        ),
    ]

    for label, points in curves:
        _print_curve(label, points)
