"""Validate independent acceleration candidates on known nonsymmetric systems."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from tests.regression.test_linear_strategy_study import capture, ROOT


class AccelerationStudyTests(unittest.TestCase):
    def test_point_directory_reset_preserves_within_point_reuse(self):
        exe=ROOT/'build-release/linear_solver_acceleration_study.exe'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'point1').mkdir();(root/'point2').mkdir()
            inputs=[root/'point1/a.bin',root/'point2/b.bin',root/'point2/c.bin']
            for i,p in enumerate(inputs):capture(p,float(i+1),True)
            out=root/'out.json'
            env=dict(os.environ,OMP_NUM_THREADS='2',VELA_LINEAR_THREADS='2',VELA_BLAS_THREADS='1',
                     VELA_STUDY_COLD_ANALYSIS='0',VELA_STUDY_RESET_PER_DIRECTORY='1')
            run=subprocess.run([str(exe),'strumpack_amd',str(out),*map(str,inputs)],env=env,capture_output=True,timeout=30)
            self.assertEqual(run.returncode,0,run.stderr.decode(errors='replace'))
            d=json.loads(out.read_text());self.assertTrue(d['reset_per_directory'])
            self.assertEqual([s['analysis_reused'] for s in d['systems']],[False,False,True])
            self.assertTrue(all(s['pass'] for s in d['systems']))

    def test_ordering_variants_cold_and_reused_analysis(self):
        exe=ROOT/'build-release/linear_solver_acceleration_study.exe'
        for mode in ('strumpack','strumpack_amd','strumpack_mmd','strumpack_and'):
            for cold in (False,True):
                with self.subTest(mode=mode,cold=cold), tempfile.TemporaryDirectory() as tmp:
                    root=Path(tmp); a=root/'a.bin'; b=root/'b.bin'; out=root/'out.json'
                    capture(a,1.,True);capture(b,2.,True)
                    env=dict(os.environ,OMP_NUM_THREADS='2',VELA_LINEAR_THREADS='2',VELA_BLAS_THREADS='1',
                             VELA_STUDY_COLD_ANALYSIS=str(int(cold)),VELA_STUDY_RESET_PER_DIRECTORY='0')
                    run=subprocess.run([str(exe),mode,str(out),str(a),str(b)],env=env,capture_output=True,timeout=30)
                    self.assertEqual(run.returncode,0,run.stderr.decode(errors='replace'))
                    d=json.loads(out.read_text());self.assertEqual(d['cold_analysis'],cold)
                    self.assertEqual(d['systems'][1]['analysis_reused'],not cold)
                    for i,row in enumerate(d['systems']):
                        self.assertTrue(row['runtime_threads_verified']);self.assertTrue(row['pass'])
                        self.assertLess(max(row['quality']['relative_solver_step_difference_by_block']),1e-10)
                        if cold or i==0:
                            self.assertGreater(row['vendor_phase_seconds']['nested dissection'],0)
                            self.assertGreater(row['vendor_phase_seconds']['symbolic factorization'],0)
                            self.assertGreaterEqual(row['reorder_other_seconds'],0)
                        else:
                            self.assertEqual(row['reorder_seconds'],0)
                            self.assertNotIn('symbolic factorization',row['vendor_phase_seconds'])

    def test_candidates_and_explicit_blas_thread_control(self):
        for suffix in ('study', 'blas_study'):
            exe = ROOT/'build-release'/f'linear_solver_acceleration_{suffix}.exe'
            self.assertTrue(exe.is_file(), 'Build both acceleration diagnostic targets')
            for mode in ('sparselu', 'umfpack', 'strumpack', 'blr6', 'blr8', 'hss6', 'hss8', 'amalg'):
                with self.subTest(suffix=suffix, mode=mode), tempfile.TemporaryDirectory() as tmp:
                    root=Path(tmp); a=root/'a.bin'; b=root/'b.bin'; out=root/'out.json'
                    capture(a,1.,True); capture(b,2.,True)
                    env=dict(os.environ, OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4',
                             VELA_LINEAR_THREADS='2', VELA_BLAS_THREADS='1',
                             VELA_STUDY_COLD_ANALYSIS='0',VELA_STUDY_RESET_PER_DIRECTORY='0')
                    run=subprocess.run([str(exe),mode,str(out),str(a),str(b)],env=env,capture_output=True,timeout=30)
                    self.assertEqual(run.returncode,0,run.stderr.decode(errors='replace'))
                    result=json.loads(out.read_text())
                    self.assertEqual(result['blas_threads'],1)
                    self.assertEqual(result['omp_threads'],1 if mode in ('sparselu','umfpack') else 2)
                    self.assertEqual(result['eigen_use_blas'],suffix=='blas_study')
                    self.assertTrue(result['systems'][1]['analysis_reused'])
                    for row in result['systems']:
                        self.assertTrue(row['pass']); self.assertNotIn('raw',row['quality'])
                        self.assertLess(max(row['quality']['relative_solver_step_difference_by_block']),1e-10)

    def test_invalid_candidate_is_rejected_and_output_preserved(self):
        exe=ROOT/'build-release/linear_solver_acceleration_study.exe'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); src=root/'a.bin'; out=root/'out.json'; capture(src,1.,True)
            args=[str(exe),'unknown',str(out),str(src)]
            self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0)
            original=out.read_bytes()
            self.assertEqual(json.loads(original)['status'],'failed')
            self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0)
            self.assertEqual(out.read_bytes(),original)


if __name__=='__main__': unittest.main()
