import importlib.util
from pathlib import Path
import unittest
from fractions import Fraction as F

spec = importlib.util.spec_from_file_location("sg_contract", Path(__file__).resolve().parents[2]/"scripts/audit_sg_representation_contract.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def samples(values, budget="1"):
    return [{"h":h,"compatibility":c,"canonical":w,"budget":budget}
            for h,c,w in values]


class RepresentationContract(unittest.TestCase):
    def test_smooth_current_values_and_increment(self):
        r=m.audit_row(samples([(-1,-2,-2),(0,0,0),(1,2,2)]),2)
        self.assertTrue(all(x["canonical_absolute_pass"] for x in r["states"]))
        self.assertTrue(r["directions"][0]["canonical_direction_pass"])

    def test_constant_gap_does_not_qualify_absolute_model(self):
        r=m.audit_row(samples([(-1,-1,9),(0,0,10),(1,1,11)]),1)
        self.assertFalse(any(x["canonical_absolute_pass"] for x in r["states"]))
        self.assertTrue(all(x["anchored_absolute_pass"] for x in r["states"]))
        self.assertTrue(r["reference_offset_is_diagnostic_only"])

    def test_rebasing_exposes_state_dependent_gap(self):
        r=m.audit_row(samples([(-1,-2,-1),(0,0,0),(1,2,1)]),1)
        self.assertEqual(sum(not x["anchored_absolute_pass"] for x in r["states"]),2)
        self.assertFalse(r["directions"][0]["endpoint_absolute_and_Jv_intervals_overlap"])

    def test_absolute_budget_can_allow_a_smooth_candidate(self):
        r=m.audit_row(samples([(-1,-2,-1),(0,0,0),(1,2,1)],"1e9"),1)
        self.assertTrue(r["directions"][0]["endpoint_absolute_and_Jv_intervals_overlap"])

    def test_sub_binary64_increment_remains_nonzero(self):
        tiny="1e-400"
        r=m.audit_row(samples([("-1e-400",0,"-1e-400"),(0,0,0),(tiny,0,tiny)],"0"),1)
        d=r["directions"][0]
        self.assertEqual(F(d["canonical_fd"]),1)
        self.assertFalse(d["compatibility_direction_pass"])
        self.assertTrue(d["nonzero_direction"])

    def test_zero_response_is_not_nonzero_evidence(self):
        r=m.audit_row(samples([(-1,0,0),(0,0,0),(1,0,0)]),0)
        self.assertFalse(r["directions"][0]["nonzero_direction"])

    def test_negative_derivative_interval_and_exact_boundary(self):
        lo,hi=m.directional_interval(-2,"0.1")
        self.assertLess(lo,hi)
        self.assertEqual(m.relative(lo,F(-2)),F("0.1"))
        self.assertEqual(m.relative(hi,F(-2)),F("0.1"))

    def test_invalid_and_incomplete_inputs_rejected(self):
        good=samples([(-1,-1,-1),(0,0,0),(1,1,1)])
        for bad in [good[:-1], good+[good[0]], samples([(-1,0,0),(0,0,0),(1,0,0)],"-1")]:
            with self.assertRaises(ValueError):m.audit_row(bad,1)
        for v in ["NaN","Infinity",float('inf'),True]:
            with self.assertRaises(ValueError):m.number(v)
        with self.assertRaises(ValueError):m.audit_row(good,1,directional_tolerance="1")


if __name__ == "__main__":unittest.main()
