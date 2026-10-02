"""Publish complete accepted CSV rows and finished VTU files incrementally."""
from __future__ import annotations
import csv
import io
import math
from pathlib import Path
import re
import shutil

from ..plt import rows_from_runner_csv, filter_rows, build_plt_document, write_plt_atomic, write_pvd_atomic
from .support import RunError


def accepted_records(path):
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except FileNotFoundError:
        return []
    # A runner may be halfway through a row. Never publish that row.
    if not text.endswith("\n"):
        text = text[:text.rfind("\n") + 1]
    return [r for r in csv.DictReader(io.StringIO(text)) if None not in r and
            all(v is not None for v in r.values()) and r.get("converged", "1") in ("1", "true")]


def vtu_metadata(path):
    data = Path(path).read_bytes()
    marker = data.find(b"<AppendedData")
    prefix = data if marker < 0 else data[:marker]
    if not data.rstrip().endswith(b"</VTKFile>"):
        return None
    text = prefix.decode("utf-8")
    result = {}
    for dtype, name, value in re.findall(r'<DataArray type="([^"]+)" Name="([^"]+)"[^>]*>\s*([^<]*)</DataArray>', text):
        if dtype == "String":
            result[name] = bytes(int(v) for v in value.split() if int(v)).decode("utf-8")
        elif value.strip():
            try:
                result[name] = float(value.strip())
            except ValueError:
                pass  # ordinary mesh arrays contain multiple values
    return result


class OutputPublisher:
    def __init__(self, root, work, result):
        self.root, self.work, self.result = Path(root), Path(work), result
        self.completed_rows = {}
        self.entries = []
        self.outputs = set()
        self.published = set()
        self.last_vtu = None
        self.last_signature = None
        self.stage_rows = []
        self.accepted_index = 0
        # Explicit FilePrefix/Time names are reserved across the entire plan.
        # Additional per-accepted-point copies must never overwrite them.
        self.reserved_names = {series["pattern"].replace("NNNN", f'{series["global_index_start"] + n:04d}')
            for stage in result.stage_plan["stages"] for series in stage["vtu_write"]["series"]
            for n in range(len(series["t"]))}

    def _path(self, name):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def update(self, index, *, final=False):
        stage = self.result.stage_plan["stages"][index]
        records = accepted_records(self.work / stage["csv"])
        if records:
            # Give WP8 a stable snapshot of complete rows. Reading the live
            # CSV twice can expose a row the runner is still writing.
            snapshot = self.work / f"stage{index:02d}_publication.csv"
            with snapshot.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
            names = [e["name"] for e in self.result.stage_plan["electrodes"]]
            sweep = stage.get("sweep", {})
            contact = sweep.get("contact", self.result.decks[index]["sweep"]["contact"])
            start, goal = stage["electrode_voltages_start_V"][contact], stage["electrode_voltages_goal_V"][contact]
            rows = rows_from_runner_csv(snapshot, sweep_contact=contact, electrodes=names,
                t_to_bias=(start, goal), electrode_voltages=stage["electrode_voltages_start_V"])
            # Wide per-contact columns are authoritative for all electrodes.
            for row, record in zip(rows, records):
                for name in names:
                    for carrier, key in (("electron", "e_current_A_per_um"), ("hole", "h_current_A_per_um"), ("total", "total_current_A_per_um")):
                        value = record.get(f"contact_{name}_{carrier}_A_per_um")
                        if value not in (None, ""):
                            row[name][key] = float(value)
                if sweep.get("mode") == "current":
                    initial = float(records[0]["current_total_A_per_um"])
                    target = sweep["stop_current_A_per_um"]
                    row["time"] = 0.0 if target == initial else (float(record["current_total_A_per_um"]) - initial) / (target - initial)
            selected = list(filter_rows(rows, stage["plt_write"]))
            self.stage_rows = rows
            path = self._path(stage["current_file"])
            prior = self.completed_rows.get(stage["current_file"], [])
            text = build_plt_document(names, prior + selected)
            signature = (index, text)
            if signature != self.last_signature:
                write_plt_atomic(path, text)
                self.outputs.add(str(path))
                self.last_signature = signature
            if final:
                self.completed_rows[stage["current_file"]] = prior + selected
        for source in sorted(self.work.glob(f"stage{index:02d}_*.vtu")):
            if str(source) in self.published:
                continue
            meta = vtu_metadata(source)
            if meta is None:
                continue
            if "PointIndex" not in meta or "t" not in meta:
                raise RunError("E602", f"Output writer failed for '{source.name}': runner metadata missing")
            t = meta["t"]
            if not math.isfinite(t):
                raise RunError("E602", f"Output writer failed for '{source.name}': non-finite t")
            visible = []
            for series in stage["vtu_write"]["series"]:
                for n, requested in enumerate(series["t"]):
                    if abs(t - requested) <= 1e-9:
                        visible.append(series["pattern"].replace("NNNN", f'{series["global_index_start"] + n:04d}'))
            if stage["vtu_write"]["per_accepted_step"]:
                prefix = Path(self.result.stage_plan["outputs"]["final_vtu"]).with_suffix("")
                prefix = str(prefix)
                if prefix.endswith("_des"):
                    prefix = prefix[:-4]
                name = f"{prefix}_{self.accepted_index:04d}_des.vtu"
                while name in self.reserved_names:
                    self.accepted_index += 1
                    name = f"{prefix}_{self.accepted_index:04d}_des.vtu"
                self.accepted_index += 1
                visible.append(name)
            if not self.result.stage_plan["outputs"]["plot_final_only"]:
                for name in dict.fromkeys(visible):
                    target = self._path(name)
                    shutil.copyfile(source, target)
                    self.outputs.add(str(target))
                    self.entries.append({"timestep": len(self.entries), "file": name})
            self.last_vtu = source
            self.published.add(str(source))
            target = self._path(self.result.stage_plan["outputs"]["final_vtu"])
            shutil.copyfile(source, target)
            self.outputs.add(str(target))
            pvd = self._path(self.result.stage_plan["outputs"]["pvd"])
            entries = self.entries or [{"timestep": 0, "file": self.result.stage_plan["outputs"]["final_vtu"]}]
            write_pvd_atomic(pvd, entries)
            self.outputs.add(str(pvd))
        if final and self.last_vtu and self.last_vtu.name.startswith(f"stage{index:02d}_"):
            for name in stage["vtu_write"]["after_stage"]:
                target = self._path(name)
                shutil.copyfile(self.last_vtu, target)
                self.outputs.add(str(target))
        return records
