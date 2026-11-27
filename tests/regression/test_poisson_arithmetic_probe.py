"""Read-only compensated-Poisson export must match the production operator."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

RUNNER = None

class PoissonArithmeticProbe(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        mesh = {
            "nodes": [{"id":i,"x":x,"y":y} for i,(x,y) in enumerate(
                [(0,0),(1,0),(0,1),(1,1),(.5,.5)])],
            "triangles": [{"id":i,"node_ids":ids,"region_id":0}
                          for i,ids in enumerate([[0,1,4],[1,3,4],[3,2,4],[2,0,4]])],
            "regions":[{"id":0,"name":"Si","material":"Si","cell_ids":[0,1,2,3]}],
            "contacts":[{"id":0,"name":"left","region_id":0,"node_ids":[0,2]},
                        {"id":1,"name":"right","region_id":0,"node_ids":[1,3]}]}
        (self.root/"mesh.json").write_text(json.dumps(mesh))
        fields=self.root/"fields";fields.mkdir()
        for key,values in {"ElectrostaticPotential":[.1,.2,.12,.18,.15],
                           "eQuasiFermiPotential":[0,0,0,0,.002],
                           "hQuasiFermiPotential":[0,0,0,0,-.003]}.items():
            with (fields/(key+"_region0.csv")).open("w",newline="") as f:
                w=csv.writer(f);w.writerow(["node_id","component0"])
                w.writerows(enumerate(values))
        self.cfg={"simulation_type":"poisson_arithmetic_model_probe",
            "mesh_file":"mesh.json","state_fields_dir":"fields",
            "output_model_json":"model.json","output_csv":"rows.csv",
            "scaling":{"mode":"unit_scaling"},
            "doping":[{"region":"Si","donors":1e17,"acceptors":0}],
            "contacts":[{"name":"left","bias":0},{"name":"right","bias":0}],
            "solver":{"method":"newton","bandgap_narrowing":"old_slotboom",
                      "poisson_residual_precision":"binary128"}}
    def tearDown(self):
        self.tmp.cleanup()
    def run_probe(self,cfg):
        p=self.root/"config.json";p.write_text(json.dumps(cfg))
        return subprocess.run([str(RUNNER),"--config",str(p)],text=True,capture_output=True,timeout=60)
    def test_export_matches_production_and_preserves_inputs(self):
        before={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in self.root.rglob("*") if p.is_file()}
        q=self.run_probe(self.cfg);self.assertEqual(q.returncode,0,q.stderr)
        status=json.loads(q.stdout);self.assertTrue(status["read_only"])
        self.assertFalse(status["qualified"]);self.assertEqual(status["accepted_updates"],0)
        model=json.loads((self.root/"model.json").read_text())
        self.assertEqual(model["schema"],"vela.poisson-arithmetic-model.v1")
        self.assertEqual(len(model["nodes"]),5)
        self.assertTrue(all(len(n)==20 for n in model["nodes"]))
        self.assertEqual(len(model["dirichlet"]),4)
        c={**self.cfg,"simulation_type":"newton_poisson_term_probe"}
        q=self.run_probe(c);self.assertEqual(q.returncode,0,q.stderr)
        with (self.root/"rows.csv").open() as f:rows=list(csv.DictReader(f))
        for row,value in zip(rows,model["zero_low_residual_double"]):
            self.assertEqual(float(row["production_residual"]),float(value))
        self.assertTrue(all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
                            for p,h in before.items()))
        original=(self.root/"model.json").read_bytes()
        q=self.run_probe(self.cfg);self.assertNotEqual(q.returncode,0)
        self.assertEqual((self.root/"model.json").read_bytes(),original)
    def test_rejects_unsupported_statistics(self):
        c=json.loads(json.dumps(self.cfg));c["solver"]["carrier_statistics"]="fermi_dirac"
        c["solver"]["poisson_residual_precision"]="double"
        q=self.run_probe(c);self.assertNotEqual(q.returncode,0)
        self.assertIn("Boltzmann",q.stderr)
        self.assertFalse((self.root/"model.json").exists())

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--runner",required=True,type=Path)
    args,rest=parser.parse_known_args();RUNNER=args.runner.resolve()
    unittest.main(argv=[sys.argv[0]]+rest)
