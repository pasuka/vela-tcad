"""Deterministic WP5 protocol fixture; never claims to solve device physics."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import time
import xml.etree.ElementTree as ET


def vtu(path, metadata, binary):
    root = ET.Element("VTKFile", type="UnstructuredGrid", version="1.0", byte_order="LittleEndian")
    grid = ET.SubElement(root, "UnstructuredGrid")
    fields = ET.SubElement(grid, "FieldData")
    for name, value in metadata.items():
        node = ET.SubElement(fields, "DataArray", type="String" if isinstance(value, str) else "Float64",
                             Name=name, NumberOfTuples="1", format="ascii")
        node.text = " ".join(str(n) for n in value.encode()) + " 0" if isinstance(value, str) else str(value)
    data = ET.tostring(root, encoding="utf-8")
    if binary:
        data = data.replace(b"</VTKFile>", b'<AppendedData encoding="raw">_\x00\x00\x00\x00</AppendedData></VTKFile>')
    path.write_bytes(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--log")
    parser.add_argument("--log-profile")
    args = parser.parse_args()
    deck = json.loads(Path(args.config).read_text())
    sweep = deck["sweep"]
    index = sweep["vtk_field_data"]["Stage"]
    mode = os.environ.get("FAKE_RUNNER_MODE", "ok")
    count = os.environ.get("FAKE_RUNNER_COUNT_FILE")
    if count:
        with Path(count).open("a") as stream:
            stream.write(json.dumps({"stage": index, "deck": deck, "pid": os.getpid(),
                "threads": os.environ.get("OMP_NUM_THREADS"), "hdf5_locking": os.environ.get("HDF5_USE_FILE_LOCKING"),
                "log_profile": args.log_profile}) + "\n")
    if mode == "reject":
        print("fixture: deck rejected", flush=True)
        return 2
    if mode == "crash":
        return 9
    start, goal = sweep["start"], sweep["stop"]
    times = [0.0] if start == goal else [0.0, .5, 1.0]
    control = sweep.get("voltage_to_current")
    if control:
        times = [0., .5, 1.]
    if mode == "partial" and index == 1:
        times = [0.0, .5]
    names = [c["name"] for c in deck["contacts"]]
    fields = ["bias_V", "converged", "current_contact", "current_total_A_per_um", "current_electron_A_per_um", "current_hole_A_per_um"]
    fields += [f"contact_{name}_{carrier}_A_per_um" for name in names for carrier in ("electron", "hole", "total")]
    with Path(deck["output_csv"]).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for point, t in enumerate(times):
            bias = start + t * (goal - start)
            if control:
                target = control["current_direction"] * control["current_points_A_per_um"][-1]
                initial = start * 1e-7
                bias = (initial + t * (target - initial)) / 1e-7
            row = {"bias_V": bias, "converged": 1, "current_contact": sweep["current_contact"],
                   "current_total_A_per_um": bias * 1e-7, "current_electron_A_per_um": bias * 1e-7, "current_hole_A_per_um": 0}
            for n, name in enumerate(names):
                current = bias * 1e-7 * (1 if n == 0 else -1)
                row.update({f"contact_{name}_electron_A_per_um": current,
                            f"contact_{name}_hole_A_per_um": 0, f"contact_{name}_total_A_per_um": current})
            writer.writerow(row); stream.flush()
            metadata = {**sweep["vtk_field_data"], "Bias_V": bias, "PointIndex": point, "t": t,
                        "GoalReached": int(t == 1 or start == goal), "goal_reached": int(t == 1 or start == goal)}
            if mode == "bad_output":
                metadata = {"Bias_V": bias}
            vtu(Path(f"{sweep['vtk_prefix']}_{point:04d}.vtu"), metadata, sweep["vtk_binary"])
            if mode == "hang":
                while True:
                    time.sleep(.1)
    Path(sweep["write_state_file"]).write_text("fixture state only\n")
    print(f"fixture stage {index} done", flush=True)
    return 1 if mode == "fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
