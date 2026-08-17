"""Regression tests for the Sentaurus ``.par`` mapping contract.

The mapping matrix is the gate that decides whether a ``.par`` parameter may
become a Vela input.  Its value comes entirely from being *strict*, so these
tests concentrate on the ways a mapping layer usually goes wrong:

* claiming coverage because a JSON field happens to exist,
* letting a parameter switch on a model the SDevice deck never selected,
* silently downgrading a model to a similar-looking one,
* and drifting out of sync with the parser and the execution IR.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import sentaurus_execution_ir as execution_ir  # noqa: E402
import sentaurus_parameter_ir as parameter_ir  # noqa: E402
import sentaurus_parameter_map as pmap  # noqa: E402

PN2D_MODELS_PAR = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "models.par"
)
BV_MODELS_PAR = (
    REPO / "reference_tcad" / "bvmethods_sentaurus2018" / "source" / "models.par"
)

# The model set the PN2D breakdown deck actually switches on.  Kept explicit
# so a test failure names the activation context rather than hiding it.
PN2D_ACTIVE_MODELS = (
    "SRH",
    "TempDependence",
    "Auger",
    "Avalanche",
    "VanOverstraeten",
    "OldSlotboom",
    "DopingDependence",
    "HighFieldSaturation",
    "eQuantumPotential",
)


class MatrixIntegrityTest(unittest.TestCase):
    """Invariants every row must satisfy, independent of any ``.par`` file."""

    def test_every_status_is_a_declared_status(self):
        for entry in pmap.MATRIX:
            self.assertIn(entry.status, pmap.STATUSES, msg=repr(entry))

    def test_every_lossy_row_explains_itself(self):
        """A non-exact status without a reason is an unreviewable claim."""
        for entry in pmap.MATRIX:
            if entry.status == pmap.STATUS_EXACT:
                continue
            self.assertTrue(
                entry.note.strip(),
                msg=f"{entry.section}.{entry.parameter} is {entry.status} "
                    "but gives no reason",
            )

    def test_rows_are_unique(self):
        """Two rows for one key would make the status depend on ordering."""
        seen: dict[tuple, pmap.MappingEntry] = {}
        for entry in pmap.MATRIX:
            self.assertNotIn(
                entry.key, seen,
                msg=f"duplicate matrix row for {entry.key}",
            )
            seen[entry.key] = entry

    def test_importable_rows_name_a_vela_target(self):
        """An importable row must say where the value lands.

        The exceptions are selectors, which the emitter consumes to choose a
        formula and never writes into a deck.  Those must say so explicitly.
        """
        for entry in pmap.MATRIX:
            if entry.status not in pmap.DEFAULT_ALLOWED_STATUSES:
                continue
            if entry.target is not None:
                continue
            self.assertTrue(
                "selector" in entry.note or "neutral" in entry.note
                or "formula selector" in entry.note,
                msg=f"{entry.section}.{entry.parameter} is importable but has "
                    "no target and is not documented as a selector",
            )


class MatrixConsistencyTest(unittest.TestCase):
    """The matrix must stay in step with the parser and the execution IR."""

    def test_supported_gates_are_known_execution_ir_models(self):
        """A gate Vela honours must be a model the execution IR can report.

        Rows whose own model is unsupported deliberately name a gate that is
        *not* in the whitelist: that is what makes them unreachable.
        """
        whitelist = set(execution_ir.SUPPORTED_MODELS)
        for entry in pmap.MATRIX:
            if entry.requires_model is None:
                continue
            if entry.status in pmap.DEFAULT_ALLOWED_STATUSES or \
                    entry.status == pmap.STATUS_APPROXIMATED:
                self.assertIn(
                    entry.requires_model, whitelist,
                    msg=f"{entry.section}.{entry.parameter} is importable but "
                        f"is gated on {entry.requires_model!r}, which the "
                        "execution IR never reports as active",
                )

    def test_parser_unsupported_sections_are_never_importable(self):
        """A section the parser calls unsupported must not sneak in here."""
        for entry in pmap.MATRIX:
            if entry.section not in parameter_ir.UNSUPPORTED_SECTIONS:
                continue
            self.assertNotIn(
                entry.status, pmap.DEFAULT_ALLOWED_STATUSES,
                msg=f"{entry.section} is declared unsupported by the parser "
                    f"but {entry.parameter} is importable in the matrix",
            )

    def test_unimplemented_ionization_models_are_never_importable(self):
        """Only selberherr and van Overstraeten-de Man exist in Vela."""
        for section in ("OkutoCrowell", "Lackner", "UniBo", "UniBo2"):
            result = pmap.classify(
                section, "a", active_models=[section],
            )
            self.assertEqual(result.status, pmap.STATUS_UNSUPPORTED_MODEL)
            self.assertFalse(result.importable_by_default)


class LookupStrictnessTest(unittest.TestCase):
    """Lookup must never let one parameter borrow another's status."""

    def test_variant_must_match_exactly(self):
        self.assertEqual(
            pmap.lookup("Bandgap", "dEg0", "OldSlotboom").status,
            pmap.STATUS_EXACT,
        )
        self.assertEqual(
            pmap.lookup("Bandgap", "dEg0", "Slotboom").status,
            pmap.STATUS_UNSUPPORTED_FORMULA,
        )

    def test_unknown_variant_does_not_match_a_variantless_row(self):
        """``Eg0(Something)`` must not inherit the plain ``Eg0`` row."""
        self.assertIsNotNone(pmap.lookup("Bandgap", "Eg0"))
        self.assertIsNone(pmap.lookup("Bandgap", "Eg0", "Something"))

    def test_formula_specific_rows_do_not_match_another_formula(self):
        self.assertEqual(
            pmap.lookup("DopingDependence", "mumin1", formula="1").status,
            pmap.STATUS_EXACT,
        )
        self.assertIsNone(
            pmap.lookup("DopingDependence", "mumin1", formula="2"),
        )

    def test_unknown_parameter_is_unmapped_not_guessed(self):
        result = pmap.classify("Bandgap", "EgSomethingNew")
        self.assertEqual(result.status, pmap.STATUS_INACTIVE)
        self.assertIsNone(result.entry)


class ActivationContextTest(unittest.TestCase):
    """A ``.par`` file is a candidate library, not an activation list."""

    def test_without_a_context_model_gated_rows_are_not_importable(self):
        result = pmap.classify("Scharfetter", "taumax", active_models=None)
        self.assertEqual(result.status, pmap.STATUS_INACTIVE)
        self.assertIn("no activated model context", result.reason)

    def test_parameters_never_enable_a_model_that_is_not_active(self):
        """The core rule: importing parameters must not turn on new physics."""
        result = pmap.classify(
            "vanOverstraetendeMan", "a", variant="low",
            active_models=["SRH"],
        )
        self.assertEqual(result.status, pmap.STATUS_INACTIVE)
        self.assertIn("must not enable a new model", result.reason)

    def test_the_same_parameter_is_importable_once_its_model_is_active(self):
        result = pmap.classify(
            "vanOverstraetendeMan", "a", variant="low",
            active_models=["SRH", "VanOverstraeten"],
        )
        self.assertEqual(result.status, pmap.STATUS_EXACT)
        self.assertTrue(result.importable_by_default)

    def test_material_constants_need_no_activation(self):
        """Permittivity is meaningful whatever the model set is."""
        result = pmap.classify("Epsilon", "epsilon", active_models=None)
        self.assertEqual(result.status, pmap.STATUS_EXACT)


class NeutralValueTest(unittest.TestCase):
    """An unrepresentable term that is switched off is not a loss."""

    def test_neutral_value_demotes_an_unsupported_row(self):
        result = pmap.classify(
            "QuantumPotentialParameters", "nu",
            active_models=["eQuantumPotential"], values=[0.0, 0.0],
        )
        self.assertEqual(result.status, pmap.STATUS_EXACT)
        self.assertIn("neutral value", result.reason)

    def test_non_neutral_value_still_blocks(self):
        result = pmap.classify(
            "QuantumPotentialParameters", "nu",
            active_models=["eQuantumPotential"], values=[0.0, 0.3],
        )
        self.assertEqual(result.status, pmap.STATUS_UNSUPPORTED_MODEL)

    def test_one_carrier_off_is_not_enough(self):
        """A term active for holes only still changes the solve."""
        result = pmap.classify(
            "HighFieldDependence", "ku",
            active_models=["HighFieldSaturation"], values=[1.0, 2.0],
        )
        self.assertEqual(result.status, pmap.STATUS_UNSUPPORTED_MODEL)

    def test_missing_values_never_demote(self):
        result = pmap.classify(
            "QuantumPotentialParameters", "nu",
            active_models=["eQuantumPotential"], values=None,
        )
        self.assertEqual(result.status, pmap.STATUS_UNSUPPORTED_MODEL)


class DoubleCountGuardTest(unittest.TestCase):
    """``dEg0`` must not be added to the gap and to ``ni`` at the same time."""

    def test_old_slotboom_dEg0_feeds_ni_only(self):
        entry = pmap.lookup("Bandgap", "dEg0", "OldSlotboom")
        self.assertEqual(entry.target, "materials[].ni")
        self.assertNotIn("bandgap_eV", entry.target)

    def test_base_gap_parameters_feed_bandgap_only(self):
        for name in ("Eg0", "alpha", "beta", "Tpar"):
            entry = pmap.lookup("Bandgap", name)
            self.assertEqual(entry.target, "materials[].bandgap_eV")

    def test_varshni_inputs_are_frozen_not_exact(self):
        """Vela stores a number, not a temperature law; say so."""
        for name in ("Eg0", "alpha", "beta", "Tpar"):
            self.assertEqual(
                pmap.lookup("Bandgap", name).status, pmap.STATUS_FROZEN,
            )


class KnownOverstatementTest(unittest.TestCase):
    """Pin the specific claims that a naive coverage count gets wrong."""

    def test_targets_are_consumed_json_keys(self):
        self.assertIn("electron_mumin1_m2_V_s", pmap.lookup(
            "DopingDependence", "mumin1", formula="1").target)
        self.assertIn("electron_field_beta", pmap.lookup(
            "HighFieldDependence", "beta0").target)
        self.assertIn("electron_saturation_velocity_m_s", pmap.lookup(
            "HighFieldDependence", "vsat0").target)
        self.assertEqual(
            "solver.electron_quantum_potential.theta",
            pmap.lookup("QuantumPotentialParameters", "theta").target,
        )

    def test_aliases_activate_canonical_mapping_gates(self):
        result = pmap.classify(
            "DopingDependence", "mumin1", formula="1",
            active_models=["DopingDep"],
        )
        self.assertEqual(pmap.STATUS_EXACT, result.status)

    def test_unknown_active_coefficient_blocks(self):
        result = pmap.classify(
            "DopingDependence", "new_coefficient", formula="3",
            active_models=["DopingDependence"],
        )
        self.assertEqual(pmap.STATUS_UNMAPPED, result.status)
        with self.assertRaises(pmap.ParameterMapError):
            pmap.assert_importable([result], allow_lossy=True)

    def test_auger_density_enhancement_blocks_the_import(self):
        """Vela's Auger is constant-coefficient; H/N0 cannot be honoured."""
        for name in ("H", "N0"):
            result = pmap.classify("Auger", name, active_models=["Auger"])
            self.assertEqual(result.status, pmap.STATUS_UNSUPPORTED_MODEL)

    def test_auger_coefficients_are_approximations_not_exact(self):
        for name in ("A", "B", "C"):
            result = pmap.classify("Auger", name, active_models=["Auger"])
            self.assertEqual(result.status, pmap.STATUS_APPROXIMATED)

    def test_constant_mobility_exponent_is_approximated(self):
        """Vela hardcodes -2.2 for both carriers; .par has 2.5 and 2.2."""
        result = pmap.classify("ConstantMobility", "Exponent")
        self.assertEqual(result.status, pmap.STATUS_APPROXIMATED)

    def test_high_field_temperature_exponents_are_approximated(self):
        for name in ("betaexp", "vsatexp"):
            result = pmap.classify(
                "HighFieldDependence", name,
                active_models=["HighFieldSaturation"],
            )
            self.assertEqual(result.status, pmap.STATUS_APPROXIMATED)

    def test_hole_quantum_potential_has_no_importable_row(self):
        """Vela has an electron quantum potential only."""
        entry = pmap.lookup("QuantumPotentialParameters", "gamma")
        self.assertIn("electron", entry.note)
        self.assertIn("electron_quantum_gamma", entry.target)


class GateBehaviourTest(unittest.TestCase):
    """``assert_importable`` must fail closed and explain itself."""

    def _results(self, *statuses):
        return [
            pmap.Classification("S", f"p{i}", None, None, status, None)
            for i, status in enumerate(statuses)
        ]

    def test_exact_and_frozen_pass(self):
        pmap.assert_importable(
            self._results(pmap.STATUS_EXACT, pmap.STATUS_FROZEN),
        )

    def test_approximated_requires_an_explicit_opt_in(self):
        results = self._results(pmap.STATUS_APPROXIMATED)
        with self.assertRaises(pmap.ParameterMapError) as ctx:
            pmap.assert_importable(results)
        self.assertIn("allow_lossy", str(ctx.exception))
        pmap.assert_importable(results, allow_lossy=True)

    def test_unsupported_is_fatal_even_with_allow_lossy(self):
        """``--allow-lossy`` accepts approximations, never missing physics."""
        for status in (pmap.STATUS_UNSUPPORTED_MODEL,
                       pmap.STATUS_UNSUPPORTED_FORMULA):
            with self.assertRaises(pmap.ParameterMapError):
                pmap.assert_importable(
                    self._results(status), allow_lossy=True,
                )

    def test_inactive_is_not_fatal(self):
        """Inactive library entries must not block an otherwise clean run."""
        pmap.assert_importable(self._results(pmap.STATUS_INACTIVE))

    def test_unmapped_is_fatal(self):
        """Unknown input in an active section must fail closed."""
        with self.assertRaises(pmap.ParameterMapError):
            pmap.assert_importable(self._results(pmap.STATUS_UNMAPPED))


class CorpusClassificationTest(unittest.TestCase):
    """End-to-end: parse a real ``.par`` and run it through the matrix."""

    @classmethod
    def setUpClass(cls):
        cls.pn2d = parameter_ir.parse_parameter_ir(PN2D_MODELS_PAR)
        cls.bv = parameter_ir.parse_parameter_ir(BV_MODELS_PAR)

    def test_every_classification_is_serialisable(self):
        results, report = pmap.classify_ir(self.pn2d, PN2D_ACTIVE_MODELS)
        self.assertTrue(results)
        json.dumps(report.to_json())

    def test_counts_add_up_to_the_parameters_examined(self):
        results, report = pmap.classify_ir(self.pn2d, PN2D_ACTIVE_MODELS)
        self.assertEqual(sum(report.counts.values()), len(results))

    def test_pn2d_blocks_on_the_auger_density_enhancement(self):
        """The one honest blocker in this corpus, named explicitly."""
        _, report = pmap.classify_ir(self.pn2d, PN2D_ACTIVE_MODELS)
        blocked = {(item.section, item.parameter) for item in report.blocking}
        self.assertEqual(blocked, {("Auger", "H"), ("Auger", "N0")})

    def test_inactive_alternative_models_do_not_block(self):
        """Lackner and UniBo ship in every ``.par`` and are never selected."""
        _, report = pmap.classify_ir(self.pn2d, PN2D_ACTIVE_MODELS)
        blocked_sections = {item.section for item in report.blocking}
        for section in ("Lackner", "UniBo", "UniBo2", "OkutoCrowell"):
            self.assertNotIn(section, blocked_sections)

    def test_bvmethods_srh_only_deck_imports_cleanly(self):
        """A file containing exactly what Vela implements must pass the gate."""
        results, report = pmap.classify_ir(
            self.bv, ("SRH", "TempDependence"),
        )
        self.assertEqual(report.blocking, [])
        self.assertEqual(report.lossy, [])
        pmap.assert_importable(results)

    def test_bvmethods_exp_temp_dependence_blocks_when_active(self):
        """``Tcoeff`` is inert under TempDep and fatal under ExpTempDep."""
        results, _ = pmap.classify_ir(
            self.bv, ("SRH", "ExpTempDependence"),
        )
        statuses = {
            item.parameter: item.status
            for item in results if item.section == "Scharfetter"
        }
        self.assertEqual(
            statuses["Tcoeff"], pmap.STATUS_UNSUPPORTED_FORMULA,
        )
        with self.assertRaises(pmap.ParameterMapError):
            pmap.assert_importable(results, allow_lossy=True)

    def test_no_activation_context_imports_nothing_model_gated(self):
        results, _ = pmap.classify_ir(self.pn2d, None)
        gated = [
            item for item in results
            if item.entry is not None and item.entry.requires_model is not None
        ]
        self.assertTrue(gated)
        for item in gated:
            self.assertEqual(item.status, pmap.STATUS_INACTIVE)

    def test_shadowed_parameters_are_not_counted(self):
        """Coverage must reflect the effective file, not every line in it."""
        ir = parameter_ir.parse_parameter_ir(PN2D_MODELS_PAR)
        blocks = ir["blocks"]
        target = next(b for b in blocks if b["section"] == "Epsilon")
        results_before, _ = pmap.classify_ir(ir, PN2D_ACTIVE_MODELS)
        target["shadowed_by"] = 9999
        results_after, _ = pmap.classify_ir(ir, PN2D_ACTIVE_MODELS)
        self.assertLess(len(results_after), len(results_before))


class MatrixSerialisationTest(unittest.TestCase):
    """The matrix must be reviewable as data, not only as code."""

    def test_matrix_round_trips_through_json(self):
        payload = pmap.matrix_as_json()
        self.assertEqual(len(payload), len(pmap.MATRIX))
        restored = json.loads(json.dumps(payload))
        self.assertEqual(restored[0]["section"], pmap.MATRIX[0].section)

    def test_every_serialised_row_carries_its_status(self):
        for row in pmap.matrix_as_json():
            self.assertIn(row["status"], pmap.STATUSES)


if __name__ == "__main__":
    unittest.main()
