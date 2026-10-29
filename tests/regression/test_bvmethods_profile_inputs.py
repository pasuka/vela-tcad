import copy
import importlib.util
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[2]/"scripts/run_bvmethods_profile_regression.py"
spec=importlib.util.spec_from_file_location("bv_profile",p)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class ProfileManifestTests(unittest.TestCase):
    def manifest(self):
        z=[0.]*2600
        return {"config":"case.json","frames":[{"label":"base","inner_V":0.,"outer_V":0.,"psi_V":z,"fn_V":z,"fp_V":z}],"outer_targets_V":[0.,.01,.1]}
    def test_physical_potentials_only(self):
        m=self.manifest();module.validate_manifest(m)
        for key in ["cached_current","primitive_fields","native_jacobian"]:
            a=copy.deepcopy(m);a[key]="oracle.json"
            with self.assertRaises(ValueError):module.validate_manifest(a)
            a=copy.deepcopy(m);a["frames"][0][key]=[1.]
            with self.assertRaises(ValueError):module.validate_manifest(a)
    def test_checkpoint_shape_and_finiteness(self):
        for bad in ["shape","nonfinite","targets"]:
            a=copy.deepcopy(self.manifest())
            if bad=="shape":a["frames"][0]["fp_V"]=[0.]
            if bad=="nonfinite":a["frames"][0]["psi_V"][0]=float("nan")
            if bad=="targets":a["outer_targets_V"]=[0.,.01,.001]
            with self.assertRaises(ValueError):module.validate_manifest(a)
    def test_incomplete_newton_run_is_not_success(self):
        m=self.manifest();s={"operator":"bvmethods_tri3_qf","frames":[{"label":"base","dimension":7801,"ledger_residual_difference_Linf":0.}],"sweep":[]}
        with self.assertRaises(ValueError):module.validate_summary(s,m)
        s["sweep"]=[{"outer_V":x,"converged":True,"inner_V":x,"current_A_per_um":0.,"equation_inf":1e-10,"load_V":0.,"trace":[{"alpha":1.},{"alpha":1.}]}for x in m["outer_targets_V"]]
        module.validate_summary(s,m)
        s["sweep"][1]["converged"]=False
        with self.assertRaises(ValueError):module.validate_summary(s,m)

if __name__=="__main__":unittest.main()
