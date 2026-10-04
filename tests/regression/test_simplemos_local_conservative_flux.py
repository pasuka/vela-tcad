import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import validate_simplemos_local_conservative_flux as flux


class DirectFeedbackTests(unittest.TestCase):
    def test_one_ampere_per_metre_uses_tcad_line_units(self):
        native=flux.native_flux_for_current(1e-6)
        particles_per_m_s=native*1e4*1e-6
        self.assertAlmostEqual(particles_per_m_s*1.602176634e-19,1.)

    def test_contact_direct_cancels_unit_feedback(self):
        r=flux.prediction(0.,1.,3.,-3.)
        self.assertEqual(r['feedback_A_per_um'],3.)
        self.assertEqual(r['total_A_per_um'],0.)

    def test_free_dipole_rejects_common_weight(self):
        self.assertEqual(flux.prediction(2.,2.,3.,0.)['total_A_per_um'],0.)

    def test_orientation_reversal_reverses_response(self):
        forward=flux.prediction(2.,3.,.01,0.)['total_A_per_um']
        reverse=flux.prediction(3.,2.,.01,0.)['total_A_per_um']
        self.assertEqual(forward,-reverse)

    def test_small_residual_after_direct_cancellation_is_not_unit_gain(self):
        r=flux.prediction(0.,.9999,1.,-1.)
        self.assertAlmostEqual(r['total_A_per_um'],-.0001)


if __name__=='__main__':unittest.main()
