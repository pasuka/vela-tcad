"""Exercise the independent strategy experiment on a known nonsymmetric system."""
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
MODES=['eigen_colamd','eigen_amd','default','unsym','noscale','pivot1',
       'given_colamd','given_amd','given_colamd_noscale','given_colamd_unsym','given_colamd_pivot1']


def capture(path,scale,solver_input=False):
    n=6
    columns=[[(j,scale*(5+j)),((j+1)%n,.3),((j+3)%n,-.2)] for j in range(n)]
    columns=[sorted(c) for c in columns]
    outer=[3*j for j in range(n+1)]
    inner=[i for c in columns for i,v in c];values=[v for c in columns for i,v in c]
    exact=[1.,-2.,3.,-4.,5.,-6.];rhs=[0.]*n
    for j,c in enumerate(columns):
        for i,v in c:rhs[i]+=v*exact[j]
    matrix=struct.pack('<QQQ',n,n,len(values))+struct.pack('<7i',*outer)+struct.pack('<18i',*inner)+struct.pack('<18d',*values)
    vector=lambda x:struct.pack('<Q',n)+struct.pack('<6d',*x)
    if solver_input:
        path.write_bytes(b'VELALU02'+struct.pack('<QQ',0,0)+matrix+vector(rhs)+vector(exact))
    else:
        path.write_bytes(b'VELALU01'+struct.pack('<QQ',0,0)+matrix+vector(rhs)+matrix+vector(rhs)+vector([1.]*n)*3+vector(exact))


class StrategyStudyTests(unittest.TestCase):
    def test_solver_input_capture_has_no_unscaled_claim(self):
        tool=ROOT/'build-release'/('linear_solver_replay'+('.exe' if os.name=='nt' else ''))
        if not tool.is_file():self.skipTest('Build linear_solver_replay')
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        env.pop('VELA_LINEAR_CAPTURE_DIR',None)
        if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;'+env.get('PATH','')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.bin';out=root/'result.json';capture(source,1.,solver_input=True)
            subprocess.run([str(tool),'sparselu',str(out),str(source)],env=env,check=True,capture_output=True)
            quality=json.loads(out.read_text())['systems'][0]['quality']
            self.assertFalse(quality['raw_available']);self.assertNotIn('raw',quality)
            self.assertEqual(quality['coordinate_scope'],'solver_input_only')
            self.assertLess(max(quality['relative_solver_step_difference_by_block']),1e-12)

    def test_all_strategies_solve_known_system_and_changed_values(self):
        name='linear_solver_strategy_study'+('.exe' if os.name=='nt' else '')
        tool=Path(os.environ.get('VELA_STRATEGY_STUDY_EXE',ROOT/'build-release'/name))
        if not tool.is_file():self.skipTest('Build the optional Release strategy diagnostic or set VELA_STRATEGY_STUDY_EXE')
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;'+env.get('PATH','')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);files=[root/'a.bin',root/'b.bin']
            for p,s in zip(files,[1.,2.]):capture(p,s)
            for mode in MODES:
                with self.subTest(mode=mode):
                    out=root/(mode+'.json')
                    subprocess.run([str(tool),mode,str(out),*map(str,files)],env=env,check=True,capture_output=True)
                    result=json.loads(out.read_text());self.assertEqual(result['status'],'pass')
                    self.assertEqual(len(result['systems']),2)
                    for row in result['systems']:
                        self.assertLess(max(row['quality']['relative_raw_step_difference_by_block']),1e-12)
                        if mode.startswith('given_') and not mode.endswith('unsym'):
                            self.assertTrue(row['matches_input_q'])


if __name__=='__main__':unittest.main()
