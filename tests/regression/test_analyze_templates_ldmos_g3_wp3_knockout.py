from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.analyze_templates_ldmos_g3_wp3_knockout import (
    assert_summary_schema_contract,
    baseline_normalized,
    canonical_hash,
    contact_gate_inputs,
    edge_signature,
    endpoint_tdrs,
    matching_physics,
    reconstructed_densities,
)
from scripts.qualify_templates_ldmos_g3_wp3_knockout import qualify
from tests.regression.test_templates_ldmos_g3_wp3_knockouts import raw_case


class AnalyzeTemplatesLdmosG3Wp3KnockoutTest(unittest.TestCase):
    def test_matching_physics_changes_exactly_registered_family(self) -> None:
        base = {
            "contacts": [{"name": "gate", "bias": 1.0}],
            "solver": {
                "mobility": {
                    "model": "constant_field",
                    "high_field_driving_force": "quasi_fermi_gradient",
                    "contact_electric_field_fallback": True,
                },
                "carrier_statistics": {"model": "fermi_dirac"},
                "bandgap_narrowing": {"model": "old_slotboom"},
            },
        }
        k1 = matching_physics(base, "k1_hfs_off")
        k2 = matching_physics(base, "k2_boltzmann")
        k3 = matching_physics(base, "k3_bgn_off")
        self.assertEqual(k1["solver"]["mobility"], {"model": "constant"})
        self.assertEqual(k1["solver"]["carrier_statistics"]["model"], "fermi_dirac")
        self.assertEqual(k2["solver"]["carrier_statistics"], {"model": "boltzmann"})
        self.assertEqual(k2["solver"]["bandgap_narrowing"]["model"], "old_slotboom")
        self.assertEqual(k3["solver"]["bandgap_narrowing"], {"model": "none"})
        self.assertEqual(k3["solver"]["carrier_statistics"]["model"], "fermi_dirac")
        self.assertEqual(base["solver"]["mobility"]["model"], "constant_field")

    def test_endpoint_tdr_selection_requires_one_exact_prefix_each(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "k1_hfs_off_vg1p0_0000_des.tdr").write_bytes(b"a")
            (root / "k1_hfs_off_vg1p166667_0000_des.tdr").write_bytes(b"b")
            (root / "k1_hfs_off_des.tdr").write_bytes(b"not an endpoint")
            result = endpoint_tdrs(root, "k1_hfs_off")
            self.assertEqual(result[1.0].name, "k1_hfs_off_vg1p0_0000_des.tdr")
            self.assertEqual(
                result[1.166666666666667].name,
                "k1_hfs_off_vg1p166667_0000_des.tdr",
            )

    def test_reconstructed_density_rejects_inconsistent_edge_copies(self) -> None:
        rows = [
            {"node0": "1", "node1": "2", "electron_density0_m3": "10",
             "electron_density1_m3": "20", "hole_density0_m3": "3",
             "hole_density1_m3": "4"},
            {"node0": "1", "node1": "3", "electron_density0_m3": "10",
             "electron_density1_m3": "30", "hole_density0_m3": "3",
             "hole_density1_m3": "5"},
        ]
        self.assertEqual(reconstructed_densities(rows)[1], (10.0, 3.0))
        rows[1]["electron_density0_m3"] = "11"
        with self.assertRaises(ValueError):
            reconstructed_densities(rows)

    def test_edge_signature_is_orientation_invariant_and_baseline_table_exact(self) -> None:
        left = [
            {"edge_id": "2", "node0": "5", "node1": "3"},
            {"edge_id": "1", "node0": "1", "node1": "2"},
        ]
        right = [
            {"edge_id": "1", "node0": "2", "node1": "1"},
            {"edge_id": "2", "node0": "3", "node1": "5"},
        ]
        self.assertEqual(edge_signature(left), edge_signature(right))
        self.assertEqual(canonical_hash([(1, 2)]), canonical_hash([(1, 2)]))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "baseline.csv"
            path.write_text(
                "bias_V,seven_node_residual_over_incident_abs_flux\n"
                "1.0,0.001\n1.166666666666667,0.002\n",
                encoding="utf-8",
            )
            self.assertEqual(baseline_normalized(path), {
                1.0: 0.001, 1.166666666666667: 0.002,
            })

    def test_contact_gate_uses_converged_boundary_state_not_bulk_sg_density(self) -> None:
        mesh = {"contacts": [{"name": "drain", "node_ids": [7]}]}
        config = {"contacts": [
            {"name": "drain", "type": "ohmic", "bias": 0.1},
            {"name": "gate", "type": "metal_gate", "bias": 1.0},
        ]}
        state = {7: {
            "psi": 0.0, "phin": 0.1, "phip": 0.1,
            "electrons_m3": 12.0, "holes_m3": 2.0,
        }}
        # Deliberately non-neutral bulk-kernel values must not contaminate the
        # independent contact boundary-state gate.
        result = contact_gate_inputs(
            mesh, config, state, {7: (99.0, 1.0)}, {7: 10.0},
        )
        self.assertEqual(result["neutrality_relative_residual"], [0.0])
        self.assertEqual(result["phin_minus_contact_V"], [0.0])
        self.assertEqual(result["phip_minus_contact_V"], [0.0])

    def test_real_qualifier_output_matches_closed_r4_schema(self) -> None:
        schema = json.loads(Path(
            "schemas/vela.templates_ldmos.g3_wp3_knockout_summary.v1.schema.json"
        ).read_text(encoding="utf-8"))
        assert_summary_schema_contract(qualify(raw_case()), schema)
        broken = qualify(raw_case())
        broken["unexpected"] = True
        with self.assertRaises(ValueError):
            assert_summary_schema_contract(broken, schema)


if __name__ == "__main__":
    unittest.main()
