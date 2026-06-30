import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class BoxFaceGeometryAuditTest(unittest.TestCase):
    def test_exports_available_edge_geometry_and_marks_missing_box_faces(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sg_edges = root / "sg_edges.csv"
            elements = root / "elements.csv"
            internal = root / "internal.csv"
            out_csv = root / "box_faces.csv"
            out_summary = root / "box_faces.md"
            dual_csv = root / "dual_faces.csv"
            dual_summary = root / "dual_faces.md"
            validation_csv = root / "true_dual_validation.csv"
            validation_summary = root / "true_dual_validation.md"

            sg_edges.write_text(
                "\n".join(
                    [
                        "point_index,bias_V,edge_id,node0,node1,x0_um,y0_um,x1_um,y1_um,edge_length_m,edge_couple_m,edge_area_proxy_m2,edge_source_integral,node0_source_integral,node1_source_integral,edge_class",
                        "0,-20,7,0,1,0,0,1,0,1e-6,5e-7,2e-12,8e-20,4e-20,4e-20,interior_bulk",
                        "0,-20,8,1,2,1,0,0,1,1.414213562373095e-6,5e-7,2e-12,0,0,0,interior_bulk",
                        "0,-20,9,2,0,0,1,0,0,1e-6,5e-7,2e-12,0,0,0,interior_bulk",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            elements.write_text(
                "\n".join(
                    [
                        "id,node0,node1,node2,region,material",
                        "3,0,1,2,R.Si,Si",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            internal.write_text(
                "\n".join(
                    [
                        "source_location_type,source_entity_id,bias_V,x_um,y_um,source_weight_or_volume_cm2_for_2D,contribution_volume_cm3_or_area_cm2_for_2D,qG_contribution_A_per_um",
                        "edge,7,-20,0.5,0,2e-8,2e-8,1e-15",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            cmd = [
                sys.executable,
                str(REPO / "scripts" / "export_box_face_geometry.py"),
                "--sg-edges",
                str(sg_edges),
                "--elements-csv",
                str(elements),
                "--internal-audit",
                str(internal),
                "--out-csv",
                str(out_csv),
                "--out-summary",
                str(out_summary),
                "--dual-out-csv",
                str(dual_csv),
                "--dual-out-summary",
                str(dual_summary),
                "--validation-out-csv",
                str(validation_csv),
                "--validation-out-summary",
                str(validation_summary),
            ]
            subprocess.run(cmd, cwd=REPO, check=True)

            with out_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            row = next(row for row in rows if row["associated_edge_id"] == "7")
            self.assertEqual(row["box_face_id"], "0")
            self.assertEqual(row["owner_node_id"], "0")
            self.assertEqual(row["neighbor_node_id"], "1")
            self.assertEqual(row["associated_edge_id"], "7")
            self.assertEqual(row["left_cell_id"], "3")
            self.assertEqual(row["right_cell_id"], "")
            self.assertAlmostEqual(float(row["edge_coupling_normal_x"]), 1.0)
            self.assertAlmostEqual(float(row["edge_coupling_normal_y"]), 0.0)
            self.assertAlmostEqual(float(row["owner_to_neighbor_direction_x"]), 1.0)
            self.assertAlmostEqual(float(row["owner_to_neighbor_direction_y"]), 0.0)
            self.assertAlmostEqual(float(row["face_normal_x"]), 1.0)
            self.assertAlmostEqual(float(row["face_normal_y"]), 0.0)
            self.assertAlmostEqual(float(row["face_tangent_x"]), -0.0)
            self.assertAlmostEqual(float(row["face_tangent_y"]), 1.0)
            self.assertAlmostEqual(float(row["face_length_cm"]), 5.0e-5)
            self.assertAlmostEqual(float(row["source_area_cm2"]), 2.0e-8)
            self.assertEqual(row["material_region"], "R.Si")
            self.assertEqual(row["boundary_type"], "interior_bulk")
            summary_text = out_summary.read_text(encoding="utf-8")
            self.assertIn("edge_coupling_normal is not a true independent dual/box face normal", summary_text)

            with dual_csv.open(newline="", encoding="utf-8") as handle:
                dual_rows = list(csv.DictReader(handle))
            self.assertEqual(len(dual_rows), 6)
            dual = dual_rows[0]
            self.assertEqual(dual["dual_face_id"], "0")
            self.assertEqual(dual["dual_type"], "median_dual")
            self.assertEqual(dual["associated_primal_edge_id"], "7")
            self.assertNotEqual(dual["dual_face_vertex0_x_um"], "")
            self.assertNotEqual(dual["dual_face_normal_x"], "")
            self.assertAlmostEqual(
                (float(dual["dual_face_normal_x"]) ** 2 + float(dual["dual_face_normal_y"]) ** 2) ** 0.5,
                1.0,
            )
            pair = [row for row in dual_rows if row["associated_primal_edge_id"] == "7"]
            self.assertEqual(len(pair), 2)
            self.assertAlmostEqual(float(pair[0]["dual_face_normal_x"]), -float(pair[1]["dual_face_normal_x"]))
            self.assertAlmostEqual(float(pair[0]["dual_face_normal_y"]), -float(pair[1]["dual_face_normal_y"]))
            self.assertIn("solver has explicit dual/sub-control-volume faces: yes", dual_summary.read_text(encoding="utf-8"))
            self.assertIn("cell_area_closure_pass: yes", validation_summary.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()





