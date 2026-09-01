import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.audit_templates_ldmos_g3_edge_inverse import (
    build_incidence,
    load_t1_contract,
    one_parameter_fit,
    permutation_p_value,
    regularization_scale,
    vector_similarity,
    weighted_minimum_norm,
)


class TemplatesLdmosG3EdgeInverseTests(unittest.TestCase):
    def test_canonical_incidence_and_cross_domain_rule(self) -> None:
        matrix = build_incidence((10, 20), [(5, 10), (10, 20), (20, 30)])
        np.testing.assert_array_equal(matrix, np.asarray([
            [-1.0, 1.0, 0.0],
            [0.0, -1.0, 1.0],
        ]))
        self.assertEqual(matrix[0, 1], -matrix[1, 1])
        self.assertEqual(np.count_nonzero(matrix[:, 0]), 1)
        self.assertEqual(np.count_nonzero(matrix[:, 2]), 1)

    def test_weighted_inverse_closes_conservation(self) -> None:
        matrix = build_incidence((10, 20), [(5, 10), (10, 20), (20, 30)])
        required = np.asarray([2.0, -3.0])
        for scale in (
            np.ones(3),
            np.asarray([0.1, 1.0, 0.2]),
            np.asarray([1.0, 0.3, 0.7]),
        ):
            result = weighted_minimum_norm(matrix, required, scale)
            np.testing.assert_allclose(matrix @ result["delta"], required, rtol=1e-12, atol=1e-12)
            self.assertLess(result["relative_closure_l2"], 1e-12)

    def test_equal_weights_are_regularization_stable(self) -> None:
        matrix = build_incidence((10, 20), [(5, 10), (10, 20), (20, 30)])
        required = np.asarray([1.0, -0.5])
        phi = np.ones(3)
        geometry = np.ones(3)
        solutions = []
        for kind in ("unweighted_l2", "abs_phi_weighted", "couple_over_length_weighted"):
            scale, _ = regularization_scale(kind, phi, geometry)
            solutions.append(weighted_minimum_norm(matrix, required, scale)["delta"])
        np.testing.assert_allclose(solutions[0], solutions[1], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(solutions[0], solutions[2], rtol=1e-12, atol=1e-12)
        self.assertGreater(vector_similarity(solutions[0], solutions[1])["cosine"], 0.999999)

    def test_one_parameter_fit_and_permutation_detect_structure(self) -> None:
        matrix = build_incidence((10, 20), [(5, 10), (6, 10), (7, 10), (20, 30), (20, 31), (20, 32)])
        phi = np.asarray([1.0, 2.0, 3.0, -1.0, -2.0, -4.0])
        labels = np.asarray([0.1, 0.4, 1.0, 0.2, 0.5, 0.9])
        required = matrix @ (0.25 * phi * labels)
        observed = one_parameter_fit(required, matrix @ (phi * labels))
        self.assertGreater(observed["explained_residual_energy"], 0.999999)
        p_value = permutation_p_value(matrix, required, phi, labels, observed["explained_residual_energy"], count=999, seed=42)
        self.assertLessEqual(p_value, 0.05)

    def test_fail_closed_on_noncanonical_or_rank_deficient_graph(self) -> None:
        with self.assertRaisesRegex(ValueError, "canonical"):
            build_incidence((10,), [(10, 5)])
        matrix = build_incidence((10, 20), [(10, 20)])
        with self.assertRaisesRegex(ValueError, "rank deficient"):
            weighted_minimum_norm(matrix, np.asarray([1.0, -1.0]), np.ones(1))

    def test_t1_contract_fails_closed_before_consuming_unknown_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "summary.json").write_text(
                json.dumps({"schema": "unknown"}), encoding="utf-8"
            )
            (root / "state_manifest.json").write_text(
                json.dumps({"schema": "unknown"}), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "summary schema mismatch"):
                load_t1_contract(root)


if __name__ == "__main__":
    unittest.main()
