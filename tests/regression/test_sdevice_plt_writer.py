"""WP8 regression tests: DF-ISE ``.plt`` and ParaView ``.pvd`` writer library.

* round trip through the legacy ``scripts/sentaurus_import.py`` reader with
  bit-identical values, and through the library's own reader;
* header block token-identical to the repository fixture layout
  (``pn2d_0v.plt`` in ``test_sentaurus_import_tools.py``) and the reader
  accepts that Sentaurus fixture;
* ``filter_rows`` for explicit t lists and ``"all_accepted_steps"``;
* ``rows_from_runner_csv`` on synthetic runner CSV files, including the
  negative sweep ``v0=0, v1=-20``, equilibrium stages, boundary-control
  columns and the terminal-balance join;
* ``.pvd`` parses with ``xml.etree`` and references existing files;
* atomic writes leave no ``.tmp`` behind and ``append_plt`` concatenates rows.
"""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import sentaurus_import  # noqa: E402
from vela_sdevice import plt  # noqa: E402

ELECTRODES = ["Cathode", "Anode"]


def electrode_values(outer, inner, e_current, h_current, total):
    return {
        "outer_voltage_V": outer,
        "inner_voltage_V": inner,
        "e_current_A_per_um": e_current,
        "h_current_A_per_um": h_current,
        "total_current_A_per_um": total,
    }


def sample_rows() -> list[dict]:
    return [
        {
            "time": 0.0,
            "Cathode": electrode_values(0.0, 0.0, -6.74115722458708e-26, 2.93597182315862e-25, 2.26185610069991e-25),
            "Anode": electrode_values(0.0, 0.0, 6.74115645827379e-26, -2.93597174652734e-25, -2.26185610069996e-25),
        },
        {
            "time": 0.1 + 0.2,
            "Cathode": electrode_values(0.0, 0.0, 1.0e-5, 2.0e16, 3.0),
            "Anode": electrode_values(0.3, 0.3, -1.0e-5, -2.0e16, -3.0),
        },
        {
            "time": 1.0,
            "Cathode": electrode_values(0.0, 0.0, 1.2345678901234567e-13, 0.0, 1.2345678901234567e-13),
            "Anode": electrode_values(1.0, 0.75, -1.2345678901234567e-13, 0.0, -1.2345678901234567e-13),
        },
    ]


# Header layout of the archived Sentaurus current file
# (tests/regression/test_sentaurus_import_tools.py, ``pn2d_0v.plt``), reduced to
# the five functions this library writes.
EXPECTED_HEADER = """
DF-ISE text

Info {
  version   = 1.0
  type      = xyplot
  datasets  = [
    "time"
    "Cathode OuterVoltage" "Cathode InnerVoltage" "Cathode eCurrent" "Cathode hCurrent" "Cathode TotalCurrent"
    "Anode OuterVoltage" "Anode InnerVoltage" "Anode eCurrent" "Anode hCurrent" "Anode TotalCurrent" ]
  functions = [
    time
    OuterVoltage InnerVoltage eCurrent hCurrent TotalCurrent
    OuterVoltage InnerVoltage eCurrent hCurrent TotalCurrent ]
}

Data {
""".lstrip()

SENTAURUS_FIXTURE = """
DF-ISE text

Info {
  version   = 1.0
  type      = xyplot
  datasets  = [
    "time"
    "Cathode OuterVoltage" "Cathode InnerVoltage" "Cathode QuasiFermiPotential" "Cathode DisplacementCurrent" "Cathode eCurrent"
    "Cathode hCurrent" "Cathode TotalCurrent" "Cathode Charge" "Anode OuterVoltage" "Anode InnerVoltage"
    "Anode QuasiFermiPotential" "Anode DisplacementCurrent" "Anode eCurrent" "Anode hCurrent" "Anode TotalCurrent"
    "Anode Charge" ]
}

Data {
      0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00  -6.74115722458708E-26
      2.93597182315862E-25   2.26185610069991E-25   5.28941943767407E-16   0.00000000000000E+00   0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   6.74115645827379E-26  -2.93597174652734E-25  -2.26185610069996E-25
     -5.28941943767407E-16
      0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00  -8.30875544197565E-26
      8.00477366113444E-25   7.17389811693687E-25   5.28941933930105E-16   0.00000000000000E+00   0.00000000000000E+00
      0.00000000000000E+00   0.00000000000000E+00   8.30875493000826E-26  -8.00477360993773E-25  -7.17389811693691E-25
     -5.28941933930104E-16
}
""".lstrip()

MAIN_CSV_COLUMNS = [
    "mode", "bias_contact", "bias_V", "current_contact", "current_total", "converged",
    "current_total_A_per_um", "current_electron_A_per_um", "current_hole_A_per_um",
]
TERMINAL_BALANCE_COLUMNS = [
    "point_index", "bias_V", "bias_contact", "contact",
    "current_electron_A_per_um", "current_hole_A_per_um", "electron_minus_hole_A_per_um",
    "current_total_A_per_um", "electron_plus_hole_A_per_um", "converged",
]


def temporary_directory(prefix: str) -> tempfile.TemporaryDirectory:
    build = REPO / "build"
    if build.is_dir():
        return tempfile.TemporaryDirectory(prefix=prefix, dir=str(build))
    return tempfile.TemporaryDirectory(prefix=prefix)


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main_csv_row(bias: float, total: float, electron=None, hole=None, *, converged: int = 1,
                 current_contact: str = "Anode") -> list[object]:
    return ["iv", "Anode", repr(bias), current_contact, repr(total * 1.0e6), converged,
            repr(total), "" if electron is None else repr(electron), "" if hole is None else repr(hole)]


def balance_row(index: int, bias: float, contact: str, electron: float, hole: float, total: float) -> list[object]:
    return [index, repr(bias), "Anode", contact, repr(electron), repr(hole), 0, repr(total), 0, 1]


class PltDocumentTest(unittest.TestCase):
    def test_round_trip_through_legacy_sentaurus_import_reader(self) -> None:
        rows = sample_rows()
        text = plt.build_plt_document(ELECTRODES, rows)

        datasets = sentaurus_import.parse_quoted_list(text, "datasets")
        self.assertEqual(datasets[0], "time")
        self.assertEqual(len(datasets), 1 + 5 * len(ELECTRODES))
        values = sentaurus_import.parse_values_block(text, len(datasets))
        self.assertEqual(len(values), len(rows))
        for row, parsed in zip(rows, values):
            expected = [row["time"]]
            for electrode in ELECTRODES:
                expected.extend(row[electrode][key] for key in plt.ROW_KEYS)
            self.assertEqual(parsed, expected)
            self.assertEqual([repr(value) for value in parsed], [repr(value) for value in expected])

        with temporary_directory("vela_plt_roundtrip_") as tmp:
            source = Path(tmp) / "pn2d_iv.plt"
            output = Path(tmp) / "pn2d_iv.csv"
            source.write_text(text, encoding="utf-8")
            sentaurus_import.export_plt_reference(source, output, "Anode OuterVoltage", "Anode TotalCurrent")
            lines = output.read_text().splitlines()
        self.assertEqual(lines[0], "bias_V,current_total")
        self.assertEqual(lines[1:], [
            f"{row['Anode']['outer_voltage_V']:.12g},{row['Anode']['total_current_A_per_um']:.12g}"
            for row in rows
        ])

    def test_round_trip_through_library_reader(self) -> None:
        rows = sample_rows()
        electrodes, parsed = plt.read_plt_document(plt.build_plt_document(ELECTRODES, rows))
        self.assertEqual(electrodes, ELECTRODES)
        self.assertEqual(parsed, rows)

    def test_header_matches_repository_fixture_layout(self) -> None:
        text = plt.build_plt_document(ELECTRODES, sample_rows())
        header = text[:text.index("Data {") + len("Data {")]
        self.assertEqual(header.split(), EXPECTED_HEADER.split())
        self.assertEqual(header.splitlines(), EXPECTED_HEADER.splitlines())
        lines = text.splitlines()
        self.assertEqual(lines[0], "DF-ISE text")
        self.assertEqual(lines[-1], "}")
        self.assertTrue(text.endswith("}\n"))
        data_lines = lines[lines.index("Data {") + 1:-1]
        self.assertEqual(len(data_lines), 3 * (1 + len(ELECTRODES)))
        self.assertEqual(data_lines[0].split(), ["0.0"])
        self.assertEqual(len(data_lines[1].split()), 5)
        self.assertEqual(data_lines[3].split(), ["0.30000000000000004"])

    def test_reader_accepts_sentaurus_fixture_with_extra_datasets(self) -> None:
        electrodes, rows = plt.read_plt_document(SENTAURUS_FIXTURE)
        self.assertEqual(electrodes, ["Cathode", "Anode"])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["time"], 0.0)
        self.assertEqual(rows[0]["Cathode"]["e_current_A_per_um"], -6.74115722458708e-26)
        self.assertEqual(rows[0]["Cathode"]["total_current_A_per_um"], 2.26185610069991e-25)
        self.assertEqual(rows[0]["Cathode"]["Charge"], 5.28941943767407e-16)
        self.assertEqual(rows[1]["Anode"]["h_current_A_per_um"], -8.00477360993773e-25)
        self.assertEqual(rows[1]["Anode"]["Charge"], -5.28941933930104e-16)
        self.assertEqual(sentaurus_import.parse_values_block(SENTAURUS_FIXTURE, 17)[1][15],
                         rows[1]["Anode"]["total_current_A_per_um"])

    def test_sign_factor_applies_to_currents_only(self) -> None:
        rows = sample_rows()
        flipped = plt.build_plt_document(ELECTRODES, rows, sign_factor=-1.0)
        _, parsed_raw = plt.read_plt_document(flipped, sign_factor=1.0)
        for row, raw in zip(rows, parsed_raw):
            self.assertEqual(raw["time"], row["time"])
            for electrode in ELECTRODES:
                self.assertEqual(raw[electrode]["outer_voltage_V"], row[electrode]["outer_voltage_V"])
                self.assertEqual(raw[electrode]["inner_voltage_V"], row[electrode]["inner_voltage_V"])
                for key in ("e_current_A_per_um", "h_current_A_per_um", "total_current_A_per_um"):
                    self.assertEqual(raw[electrode][key], -row[electrode][key])
        _, restored = plt.read_plt_document(flipped, sign_factor=-1.0)
        self.assertEqual(restored, rows)
        self.assertEqual(plt.CURRENT_SIGN_FACTOR, 1.0)
        self.assertEqual(plt.build_plt_document(ELECTRODES, rows),
                         plt.build_plt_document(ELECTRODES, rows, sign_factor=plt.CURRENT_SIGN_FACTOR))

    def test_none_currents_are_written_as_zero(self) -> None:
        row = {"time": 0.5, "Anode": electrode_values(1.0, 1.0, None, None, 2.0e-9)}
        _, parsed = plt.read_plt_document(plt.build_plt_document(["Anode"], [row]))
        self.assertEqual(parsed[0]["Anode"]["e_current_A_per_um"], 0.0)
        self.assertEqual(parsed[0]["Anode"]["h_current_A_per_um"], 0.0)
        self.assertEqual(parsed[0]["Anode"]["total_current_A_per_um"], 2.0e-9)

    def test_invalid_input_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            plt.build_plt_document(["Anode", "Anode"], [])
        with self.assertRaises(ValueError):
            plt.build_plt_document(["bad name"], [])
        with self.assertRaises(ValueError):
            plt.build_plt_document(["Anode"], [{"time": float("nan"), "Anode": electrode_values(0, 0, 0, 0, 0)}])
        with self.assertRaises(ValueError):
            plt.build_plt_document(["Anode"], [{"time": 0.0, "Anode": {"outer_voltage_V": 0.0}}])
        with self.assertRaises(ValueError):
            plt.read_plt_document("Info {\n datasets = [ \"time\" \"a OuterVoltage\" ]\n}\nData {\n 0 1 2\n}\n")


class FilterRowsTest(unittest.TestCase):
    def test_filter_rows_by_t_list_and_all_accepted_steps(self) -> None:
        rows = [{"time": t} for t in (0.0, 0.25, 0.5, 0.75, 1.0)]
        self.assertEqual(plt.filter_rows(rows, "all_accepted_steps"), rows)
        self.assertEqual(plt.filter_rows(rows, [0.0, 0.5 + 5.0e-10, 1.0]), [rows[0], rows[2], rows[4]])
        self.assertEqual(plt.filter_rows(rows, [0.5 + 2.0e-9]), [])
        self.assertEqual(plt.filter_rows(rows, []), [])
        with self.assertRaises(ValueError):
            plt.filter_rows(rows, "every_step")


class RowsFromRunnerCsvTest(unittest.TestCase):
    def test_negative_sweep_with_terminal_balance(self) -> None:
        with temporary_directory("vela_plt_runner_csv_") as tmp:
            root = Path(tmp)
            main = root / "stage01.csv"
            balance = root / "stage01_terminal_balance.csv"
            write_csv(main, MAIN_CSV_COLUMNS, [
                main_csv_row(0.0, -1.0e-17, -1.5e-17, 0.5e-17),
                main_csv_row(-5.0, -2.0e-16, -2.5e-16, 0.5e-16),
                main_csv_row(-7.5, -3.0e-16, -3.5e-16, 0.5e-16, converged=0),
                main_csv_row(-10.0, -4.0e-16),
                main_csv_row(-20.0, -8.0e-16, -8.5e-16, 0.5e-16),
            ])
            write_csv(balance, TERMINAL_BALANCE_COLUMNS, [
                balance_row(0, 0.0, "Anode", -1.5e-17, 0.5e-17, -1.0e-17),
                balance_row(0, 0.0, "Cathode", 1.5e-17, -0.5e-17, 1.0e-17),
                balance_row(1, -5.0, "Anode", -2.5e-16, 0.5e-16, -2.0e-16),
                balance_row(1, -5.0, "Cathode", 2.5e-16, -0.5e-16, 2.0e-16),
                balance_row(2, -10.0, "Anode", -4.5e-16, 0.5e-16, -4.0e-16),
                balance_row(2, -10.0, "Cathode", 4.5e-16, -0.5e-16, 4.0e-16),
                balance_row(3, -20.0, "Anode", -8.5e-16, 0.5e-16, -8.0e-16),
                balance_row(3, -20.0, "Cathode", 8.5e-16, -0.5e-16, 8.0e-16),
            ])
            rows = plt.rows_from_runner_csv(
                main,
                sweep_contact="Anode",
                electrodes=ELECTRODES,
                t_to_bias=(0.0, -20.0),
                terminal_balance_csv_path=balance,
                electrode_voltages={"Cathode": 0.0, "Anode": 0.0},
            )

        self.assertEqual([row["time"] for row in rows], [0.0, 0.25, 0.5, 1.0])
        self.assertEqual([repr(row["time"]) for row in rows], ["0.0", "0.25", "0.5", "1.0"])
        self.assertEqual([row["Anode"]["outer_voltage_V"] for row in rows], [0.0, -5.0, -10.0, -20.0])
        self.assertEqual([row["Anode"]["inner_voltage_V"] for row in rows], [0.0, -5.0, -10.0, -20.0])
        self.assertEqual([row["Cathode"]["outer_voltage_V"] for row in rows], [0.0] * 4)
        self.assertEqual([row["Anode"]["total_current_A_per_um"] for row in rows],
                         [-1.0e-17, -2.0e-16, -4.0e-16, -8.0e-16])
        self.assertEqual(rows[1]["Anode"]["e_current_A_per_um"], -2.5e-16)
        self.assertEqual(rows[1]["Anode"]["h_current_A_per_um"], 0.5e-16)
        # Carrier split missing from the main CSV falls back to the terminal balance.
        self.assertEqual(rows[2]["Anode"]["e_current_A_per_um"], -4.5e-16)
        self.assertEqual(rows[2]["Anode"]["h_current_A_per_um"], 0.5e-16)
        self.assertEqual([row["Cathode"]["total_current_A_per_um"] for row in rows],
                         [1.0e-17, 2.0e-16, 4.0e-16, 8.0e-16])
        self.assertEqual(rows[3]["Cathode"]["e_current_A_per_um"], 8.5e-16)
        self.assertEqual(rows[3]["Cathode"]["h_current_A_per_um"], -0.5e-16)
        self.assertEqual(plt.filter_rows(rows, [0.5, 1.0]), rows[2:])
        plt.read_plt_document(plt.build_plt_document(ELECTRODES, rows))

    def test_main_csv_only_leaves_none_carriers_and_zero_other_electrodes(self) -> None:
        with temporary_directory("vela_plt_runner_csv_main_") as tmp:
            main = Path(tmp) / "stage01.csv"
            write_csv(main, ["bias_V", "bias_contact", "current_contact", "current_total_A_per_um", "converged"], [
                [repr(0.0), "Anode", "Anode", repr(1.0e-18), 1],
                [repr(0.3), "Anode", "Anode", repr(2.0e-12), 1],
                [repr(0.6), "Anode", "Anode", repr(3.0e-6), 1],
            ])
            rows = plt.rows_from_runner_csv(
                main, sweep_contact="Anode", electrodes=ELECTRODES, t_to_bias=(0.0, 0.6),
                electrode_voltages={"Cathode": 0.0})
        self.assertEqual([row["time"] for row in rows], [0.0, 0.5, 1.0])
        self.assertIsNone(rows[1]["Anode"]["e_current_A_per_um"])
        self.assertIsNone(rows[1]["Anode"]["h_current_A_per_um"])
        self.assertEqual(rows[1]["Anode"]["total_current_A_per_um"], 2.0e-12)
        self.assertEqual(rows[1]["Cathode"], electrode_values(0.0, 0.0, 0.0, 0.0, 0.0))
        _, parsed = plt.read_plt_document(plt.build_plt_document(ELECTRODES, rows))
        self.assertEqual(parsed[1]["Anode"]["e_current_A_per_um"], 0.0)
        self.assertEqual(parsed[1]["Anode"]["total_current_A_per_um"], 2.0e-12)

    def test_equilibrium_stage_has_time_zero_and_fixed_voltages(self) -> None:
        with temporary_directory("vela_plt_runner_csv_eq_") as tmp:
            main = Path(tmp) / "stage00.csv"
            write_csv(main, MAIN_CSV_COLUMNS, [main_csv_row(0.7, 1.0e-9, 1.5e-9, -0.5e-9)])
            rows = plt.rows_from_runner_csv(
                main, sweep_contact="Anode", electrodes=["Anode", "Base", "Cathode"], t_to_bias=(0.7, 0.7),
                electrode_voltages={"Base": 0.7, "Cathode": 3.0})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["time"], 0.0)
        self.assertEqual(rows[0]["Anode"]["outer_voltage_V"], 0.7)
        self.assertEqual(rows[0]["Base"]["outer_voltage_V"], 0.7)
        self.assertEqual(rows[0]["Cathode"]["inner_voltage_V"], 3.0)
        self.assertEqual(rows[0]["Anode"]["e_current_A_per_um"], 1.5e-9)

    def test_boundary_control_columns_and_column_map(self) -> None:
        with temporary_directory("vela_plt_runner_csv_bc_") as tmp:
            main = Path(tmp) / "stage02.csv"
            write_csv(main, ["bias_V", "current_contact", "I_total", "ok", "inner_voltage_V", "outer_voltage_V"], [
                [repr(5.9), "drain", repr(1.0e-7), "1", repr(5.9), repr(6.9)],
                [repr(6.0), "drain", repr(2.0e-7), "0", repr(6.0), repr(8.0)],
                [repr(6.1), "drain", repr(3.0e-7), "1", repr(6.1), repr(9.1)],
            ])
            rows = plt.rows_from_runner_csv(
                main, sweep_contact="drain", electrodes=["drain", "source"], t_to_bias=(6.9, 9.1),
                column_map={"total_current_A_per_um": "I_total", "converged": "ok"})
            with self.assertRaises(ValueError):
                plt.rows_from_runner_csv(main, sweep_contact="drain", electrodes=["drain", "source"],
                                         t_to_bias=(6.9, 9.1))
            with self.assertRaises(ValueError):
                plt.rows_from_runner_csv(main, sweep_contact="gate", electrodes=["drain", "source"],
                                         t_to_bias=(6.9, 9.1), column_map={"total_current_A_per_um": "I_total"})
        self.assertEqual([row["time"] for row in rows], [0.0, 1.0])
        self.assertEqual([row["drain"]["outer_voltage_V"] for row in rows], [6.9, 9.1])
        self.assertEqual([row["drain"]["inner_voltage_V"] for row in rows], [5.9, 6.1])
        self.assertEqual([row["drain"]["total_current_A_per_um"] for row in rows], [1.0e-7, 3.0e-7])


class FileWritersTest(unittest.TestCase):
    def test_pvd_parses_and_references_existing_files(self) -> None:
        with temporary_directory("vela_plt_pvd_") as tmp:
            root = Path(tmp)
            entries = []
            for index in (0, 1, 2, 5):
                name = f"pn2d_iv_multibias_{index:04d}_des.vtu"
                (root / name).write_text("", encoding="utf-8")
                entries.append({"timestep": index, "file": name})
            pvd = root / "pn2d_iv_multibias.pvd"
            plt.write_pvd_atomic(pvd, entries)
            self.assertFalse((root / (pvd.name + ".tmp")).exists())
            self.assertTrue(pvd.read_text(encoding="utf-8").startswith('<?xml version="1.0"?>\n'))

            tree = ET.parse(pvd)
            root_element = tree.getroot()
            self.assertEqual(root_element.tag, "VTKFile")
            self.assertEqual(root_element.attrib["type"], "Collection")
            datasets = root_element.findall("./Collection/DataSet")
            self.assertEqual([int(node.attrib["timestep"]) for node in datasets], [0, 1, 2, 5])
            for node in datasets:
                self.assertEqual(node.attrib["group"], "")
                self.assertEqual(node.attrib["part"], "0")
                self.assertTrue((pvd.parent / node.attrib["file"]).is_file(), node.attrib["file"])
        with self.assertRaises(ValueError):
            plt.build_pvd_document([{"timestep": 0, "file": "a.vtu"}, {"timestep": 0, "file": "b.vtu"}])
        with self.assertRaises(ValueError):
            plt.build_pvd_document([{"timestep": "0", "file": "a.vtu"}])
        with self.assertRaises(ValueError):
            plt.build_pvd_document([{"timestep": 1, "file": ""}])

    def test_atomic_write_and_append(self) -> None:
        rows = sample_rows()
        with temporary_directory("vela_plt_atomic_") as tmp:
            root = Path(tmp)
            target = root / "pn2d_iv.plt"
            text = plt.build_plt_document(ELECTRODES, rows[:1])
            plt.write_plt_atomic(target, text)
            self.assertEqual(target.read_bytes(), text.encode("utf-8"))
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["pn2d_iv.plt"])

            written = plt.append_plt(target, ELECTRODES, rows[1:2])
            self.assertEqual(written, rows[:2])
            written = plt.append_plt(target, ELECTRODES, rows[2:])
            self.assertEqual(written, rows)
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["pn2d_iv.plt"])
            self.assertEqual(target.read_bytes(), plt.build_plt_document(ELECTRODES, rows).encode("utf-8"))
            electrodes, parsed = plt.read_plt_document(target.read_text(encoding="utf-8"))
            self.assertEqual(electrodes, ELECTRODES)
            self.assertEqual(parsed, rows)

            with self.assertRaises(ValueError):
                plt.append_plt(target, ["Anode", "Cathode"], rows[:1])
            self.assertEqual(target.read_bytes(), plt.build_plt_document(ELECTRODES, rows).encode("utf-8"))

            fresh = root / "fresh_des.plt"
            plt.append_plt(fresh, ["Anode"], [{"time": 0.0, "Anode": electrode_values(0, 0, 0, 0, 0)}])
            self.assertEqual(plt.read_plt_document(fresh.read_text(encoding="utf-8"))[0], ["Anode"])
            self.assertFalse((root / "fresh_des.plt.tmp").exists())

    def test_failed_write_leaves_previous_file_and_no_tmp(self) -> None:
        with temporary_directory("vela_plt_atomic_fail_") as tmp:
            root = Path(tmp)
            target = root / "pn2d_iv.plt"
            plt.write_plt_atomic(target, "original\n")
            with self.assertRaises(TypeError):
                plt.write_plt_atomic(target, None)  # type: ignore[arg-type]
            self.assertEqual(target.read_text(encoding="utf-8"), "original\n")
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["pn2d_iv.plt"])

            plt.write_plt_atomic(target, plt.build_plt_document(ELECTRODES, sample_rows()[:1]))
            with self.assertRaises(ValueError):
                plt.append_plt(target, ELECTRODES,
                               [{"time": float("inf"), "Cathode": electrode_values(0, 0, 0, 0, 0),
                                 "Anode": electrode_values(0, 0, 0, 0, 0)}])
            self.assertEqual(plt.read_plt_document(target.read_text(encoding="utf-8"))[1], sample_rows()[:1])
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["pn2d_iv.plt"])


class VerifiedCurrentSignTests(unittest.TestCase):
    def test_first_three_verified_native_anode_points_freeze_global_sign(self):
        # T-2022.03-SP2 v2 pn2d_iv.plt SHA-256:
        # a5f712e3ba0bfeea221e0de681ce380b3ad4ab14eb42c63b1324eac431709b25
        # Raw samples are log-current interpolation of accepted Anode CSV rows;
        # they establish direction only, not magnitude agreement at these biases.
        from vela_sdevice.plt.document import CURRENT_SIGN_FACTOR
        evidence = [(0.01,5.80543692864973e-19,5.3887968315032304e-18),
                    (0.02284,1.62155570626151e-18,1.3977546608561714e-17),
                    (0.03942928,3.74720553099733e-18,2.9009156046579526e-17)]
        self.assertEqual(CURRENT_SIGN_FACTOR,1.0)
        for bias,native,raw in evidence:
            text=plt.build_plt_document(['Anode'],[{'time':bias/10.,'Anode':electrode_values(bias,bias,raw,0.,raw)}])
            _,rows=plt.read_plt_document(text,sign_factor=1.)
            self.assertGreater(rows[0]['Anode']['total_current_A_per_um']*native,0.)
            self.assertEqual(rows[0]['Anode']['total_current_A_per_um'],raw*CURRENT_SIGN_FACTOR)


if __name__ == "__main__":
    unittest.main()
