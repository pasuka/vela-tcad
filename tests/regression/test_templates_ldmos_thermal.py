"""Thermal acceptance gates must not pass partial or non-finite evidence."""
import json
import math
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_templates_ldmos_thermal import assess


class ThermalGates(unittest.TestCase):
    def setUp(self):
        self.contract=json.loads((ROOT/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json').read_text(encoding='utf-8'))

    def test_absolute_floor_and_balance(self):
        r=assess([300.,301.],[300.5,301.5],[1.,3.],100.,100.05,self.contract)
        self.assertEqual(r['status'],'pass')
        self.assertEqual(r['gates']['temperature_rise_field']['rms_error_K'],.5)
        self.assertEqual(r['gates']['temperature_rise_field']['limit_K'],1.)
        self.assertEqual(assess([300.],[300.],[1.],100.,101.,self.contract)['status'],'fail')

    def test_field_gate_detects_wrong_spatial_distribution_with_equal_peak(self):
        r=assess([400.,300.],[300.,400.],[1.,3.],1.,1.,self.contract)
        self.assertEqual(r['gates']['peak_temperature_rise']['status'],'pass')
        self.assertEqual(r['gates']['temperature_rise_field']['rms_error_K'],100.)
        self.assertEqual(r['gates']['temperature_rise_field']['reference_rms_rise_K'],50.)
        self.assertEqual(r['status'],'fail')

    def test_relative_peak_threshold(self):
        self.assertEqual(assess([400.],[405.],[1.],1.,1.,self.contract)['status'],'pass')
        self.assertEqual(assess([400.],[405.01],[1.],1.,1.,self.contract)['status'],'fail')

    def test_zero_power_is_not_relative_acceptance(self):
        self.assertEqual(assess([300.],[300.],[1.],0.,0.,self.contract)['status'],'not_evaluable')

    def test_incomplete_or_nonfinite_evidence_is_rejected(self):
        for reference,candidate,volume in [([],[],[]),([300.],[300.,301.],[1.]),([300.],[math.nan],[1.]),([300.],[300.],[0.]),([300.],[300.],[-1.])]:
            with self.assertRaises(ValueError):assess(reference,candidate,volume,1.,1.,self.contract)
        with self.assertRaises(ValueError):assess([300.],[300.],[1.],math.nan,1.,self.contract)


if __name__=='__main__':unittest.main()
