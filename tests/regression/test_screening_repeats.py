import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"scripts"))
from run_templates_ldmos_screening_repeats import METHODS, FIELDS, select_candidate, paired_order, compare_states


class ScreeningRepeatProtocolTest(unittest.TestCase):
    def batch(self):
        return dict(status="completed", points=31,
                    runs=[dict(method=m, gate=g, wall_seconds=100-i, original_point_gates_pass=True)
                          for i, m in enumerate(METHODS) for g in (4, 8)],
                    comparisons=[dict(method=m, gate=g, point=p, passed=True)
                                 for m in METHODS[1:] for g in (4, 8) for p in range(31)])

    def test_fastest_must_pass_both_native_and_every_state_check(self):
        batch=self.batch(); native=dict.fromkeys(METHODS,"pass")
        self.assertEqual(select_candidate(batch,native),"toms748")
        native["toms748"]="fail"
        self.assertEqual(select_candidate(batch,native),"halley")
        next(c for c in batch["comparisons"] if c["method"]=="halley")["passed"]=False
        self.assertEqual(select_candidate(batch,native),"newton")
        native["legacy"]="fail"
        with self.assertRaises(ValueError):select_candidate(batch,native)

    def test_duplicate_or_missing_evidence_cannot_qualify(self):
        batch=self.batch();batch["runs"][-1]=copy.deepcopy(batch["runs"][-2])
        with self.assertRaises(ValueError):select_candidate(batch,dict.fromkeys(METHODS,"pass"))
        batch=self.batch();batch["comparisons"]=[c for c in batch["comparisons"] if c["point"]!=30]
        with self.assertRaises(ValueError):select_candidate(batch,dict.fromkeys(METHODS,"pass"))

    def test_order_balances_first_position_and_alternates_rounds(self):
        orders=[paired_order(r,g,"halley") for r in range(3) for g in (4,8)]
        self.assertEqual(sum(o[0]=="halley" for o in orders),3)
        for g in (4,8):self.assertEqual(paired_order(0,g,"halley"),tuple(reversed(paired_order(1,g,"halley"))))

    def test_original_state_gate_rejects_nonfinite_and_out_of_gate_values(self):
        state={k:[0.,1.] for k in FIELDS}
        self.assertEqual(max(compare_states(state,state).values()),0.)
        for invalid in (float("nan"),float("inf"),2e-8):
            changed=copy.deepcopy(state);changed[FIELDS[0]][0]=invalid
            with self.assertRaises(ValueError):compare_states(state,changed)


if __name__=="__main__":unittest.main()
