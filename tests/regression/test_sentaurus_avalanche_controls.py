#!/usr/bin/env python3
"""Tests for explicit Sentaurus avalanche driving-force control decks."""

from __future__ import annotations

import unittest
from unittest import mock

from scripts.sentaurus_avalanche_controls import (
    VARIANTS,
    make_variant_deck,
    sentaurus_release,
    validate_remote_root,
)


BASE_DECK = """File {
  Plot = "runtime_element_avalanche_probe_default.tdr"
  Current = "runtime_element_avalanche_probe_default.plt"
  Output = "runtime_element_avalanche_probe_default"
}
Physics {
  Mobility(DopingDependence HighFieldSaturation)
  Recombination(SRH Avalanche(VanOverstraeten))
}
Math {
  Extrapolate
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""


class SentaurusAvalancheDriveControlsTest(unittest.TestCase):
    def test_variant_contract_covers_four_orthogonal_branches(self) -> None:
        self.assertEqual(
            tuple(VARIANTS),
            (
                "implicit_default",
                "explicit_grad_qf",
                "explicit_electric_field",
                "grad_qf_aval_dens_grad_qf",
            ),
        )

    def test_decks_select_only_the_declared_avalanche_controls(self) -> None:
        decks = {
            name: make_variant_deck(BASE_DECK, name, (-1, -10, -20))
            for name in VARIANTS
        }
        self.assertIn(
            "Avalanche(VanOverstraeten)",
            decks["implicit_default"],
        )
        self.assertNotIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            decks["implicit_default"],
        )
        self.assertNotIn("AvalDensGradQF", decks["implicit_default"])

        self.assertIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            decks["explicit_grad_qf"],
        )
        self.assertNotIn("AvalDensGradQF", decks["explicit_grad_qf"])

        self.assertIn(
            "Avalanche(VanOverstraeten ElectricField)",
            decks["explicit_electric_field"],
        )
        self.assertNotIn(
            "GradQuasiFermi",
            decks["explicit_electric_field"],
        )

        combined = decks["grad_qf_aval_dens_grad_qf"]
        self.assertIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            combined,
        )
        self.assertEqual(combined.count("AvalDensGradQF"), 1)

    def test_decks_have_unique_outputs_and_exact_bias_goals(self) -> None:
        for name in VARIANTS:
            deck = make_variant_deck(BASE_DECK, name, (-1, -10, -20))
            self.assertIn(
                f'runtime_element_avalanche_probe_{name}.plt',
                deck,
            )
            for bias in (-1, -10, -20):
                self.assertEqual(
                    deck.count(f'Goal {{ Name="Anode" Voltage={bias} }}'),
                    1,
                )

    def test_remote_root_rejects_shell_syntax_and_non_normal_paths(self) -> None:
        self.assertEqual(
            validate_remote_root("/home/tcad/safe_root"),
            "/home/tcad/safe_root",
        )
        for invalid in (
            "relative/path",
            "/home/tcad/has space",
            "/home/tcad/root;touch_bad",
            "/home/tcad/../bad",
            "/home//tcad",
            "/home/tcad/",
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    validate_remote_root(invalid)

    def test_release_probe_accepts_sdevice_nonzero_version_exit(self) -> None:
        completed = mock.Mock(
            returncode=1,
            stdout=(
                "path=/opt/sentaurus/sdevice\n"
                "*** Version O-2018.06-SP2 ***\n"
            ),
            stderr="",
        )
        with mock.patch("subprocess.run", return_value=completed):
            self.assertEqual(
                sentaurus_release("ssh", "sentaurus"),
                "O-2018.06-SP2",
            )


if __name__ == "__main__":
    unittest.main()
