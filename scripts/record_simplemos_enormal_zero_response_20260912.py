"""Keep the nonresponsive native stress-parameter experiment explicitly failed."""
from pathlib import Path
import check_simplemos_enormal_response_20260912 as q
a,d,L,O=q.a,q.d,q.L,q.O
def main():
    a.verify(O/'zero_evidence.json');rows=a.rows(O/'signed_points.csv');out=[]
    for case in sorted({r['case'] for r in rows}):
        group=[r for r in rows if r['case']==case];ref=float(next(r for r in group if r['label']=='zero')['Id_A_per_um'])
        assert all(r['native_qualified']=='True' for r in group)
        change=max(abs(float(r['Id_A_per_um'])/ref-1) for r in group)
        assert change==0
        out.append(dict(case=case,points=len(group),max_current_change_relative=change,response_qualified=False,reason='a_ac and a_sr are stress enhancement factors, inactive without the corresponding stress model; zero signal cannot calibrate direct inverse-scattering scaling.'))
    a.write_csv(O/'failed_response.csv',out)
    a.write(O/'failure_summary.json',dict(native_exit_and_bias_qualified=20,response_qualified=0,acceptance_relaxed=False,analysis_error='Original response analyzer divided by the exactly zero central slope; preserve its traceback. This record classifies the zero signal as failure, not a linearity pass.',manual='T-2022.03 Device User Guide page 995, Factor Models Applied to Mobility Components, Equation 1074.'))
    d.matrix.freeze(O/'failure_evidence.json',[Path(__file__).resolve(),O/'zero_evidence.json',O/'signed_points.csv',O/'failed_response.csv',O/'failure_summary.json',L/'signed_results.tgz',L/'signed_check.log']+[f for f in (L/'signed_raw').rglob('*') if f.is_file()])
    print(out,flush=True)
if __name__=='__main__':main()
