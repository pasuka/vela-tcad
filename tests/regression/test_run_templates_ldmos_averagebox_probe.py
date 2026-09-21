import unittest

from scripts.run_templates_ldmos_averagebox_probe import (
    render_deck,
    render_state_capture_deck,
)


class TemplatesLdmosAverageBoxProbeTest(unittest.TestCase):
    def test_render_deck_adds_both_debug_contracts_and_one_step_solve(self) -> None:
        source = """Plot {\n  eDensity\n}\nMath {\n  Iterations=25\n}\nSolve {\n  Coupled { Poisson Electron Hole }\n}\n"""
        rendered = render_deck(source)
        self.assertIn("AverageBoxMethod", rendered)
        self.assertIn("BoxMeasureFromFile(GrdNumbering)", rendered)
        self.assertIn("BoxCoefficientsFromFile(GrdNumbering)", rendered)
        self.assertIn("BM_CoeffIntersectionNonDelaunayElements", rendered)
        self.assertIn("Coupled(Iterations=1) { Poisson }", rendered)
        self.assertNotIn("Poisson Electron Hole", rendered)

    def test_render_deck_rejects_preexisting_box_policy(self) -> None:
        with self.assertRaises(ValueError):
            render_deck("Plot {\n}\nMath { AverageBoxMethod\n}\nSolve {\n}\n")

    def test_render_deck_can_load_state_for_one_coupled_iteration(self) -> None:
        source = """Electrode {
  { name= "gate" Material= "PolySi"(N) Voltage= 0.0}
}
Plot {
  eDensity
}
Math {
  Iterations=25
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""
        rendered = render_deck(
            source, initial_state_prefix="initial_state", gate_voltage_V=0.5
        )
        self.assertIn('Voltage= 0.5', rendered)
        self.assertIn('Load(FilePrefix="initial_state")', rendered)
        self.assertIn("Coupled(Iterations=1) { Poisson Electron Hole }", rendered)

    def test_render_state_capture_deck_changes_goal_and_saves(self) -> None:
        source = """Solve {
  Quasistationary(Goal { name= "gate" voltage= 1.0 }) {
    Coupled { Poisson Electron Hole }
  }
}
"""
        rendered = render_state_capture_deck(source, 0.5)
        self.assertIn('Goal { name= "gate" voltage= 0.5 }', rendered)
        self.assertIn('Save(FilePrefix="initial_state")', rendered)
        self.assertLess(rendered.index("Save("), rendered.rfind("}"))


if __name__ == "__main__":
    unittest.main()
