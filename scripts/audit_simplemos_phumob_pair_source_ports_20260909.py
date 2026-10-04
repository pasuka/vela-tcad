"""Independent current-gradient/adjoint validation of the frozen pair source."""
import math
from pathlib import Path
import validate_simplemos_phumob_pair_source_scaled_20260909 as checked

s=checked.s
OUT=s.OUT/'port_calibration'

def run():
    s.a.verify(s.OUT/'evidence.json');scope=s.a.rows(s.OUT/'scope.csv');tangent=s.a.rows(s.OUT/'tangent.csv');ports=s.a.rows(s.OUT/'ports.csv');files=[Path(__file__).resolve(),s.OUT/'evidence.json']
    cases=s.a.read(s.OUT/'contract.json')['cases']
    for c in cases:files.append(s.LOCAL/'fixed'/c['key']/'zero/adjoint.json')
    s.a.write(OUT/'contract.json',dict(scope='Same frozen pair source, direct physical current-gradient dot full-J tangent, adjoint dot source, and actual signed DC Id.',gates=dict(prediction=1e-3,duality=1e-6,two_amplitude=1e-3,even_over_odd=.01,signal_over_drift=100),does_not_override_field_response_failures=True))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);rows=[]
    for c in cases:
        path=s.LOCAL/'fixed'/c['key']/'zero/adjoint.json';status=s.execute(path,c,0.);assert status['exit_code']==0,status
        adj=s.a.rows(path.with_suffix('.csv'));tr=[r for r in tangent if r['key']==c['key']]
        direct=math.fsum(float(adj[int(r['node'])]['dI_d'+k+'_A_per_um_per_V'])*float(r[k+'_V']) for r in tr for k in ('psi','phin','phip'))
        dual=-math.fsum((float(adj[int(r['node'])]['lambda_electron'])+float(adj[int(r['node'])]['lambda_hole']))*float(r['internal_source']) for r in scope if r['key']==c['key'])
        dual_error=abs(direct-dual)/max(abs(direct),abs(dual),1e-300)
        root=s.LOCAL/'dc'/c['key'];curr=lambda label:s.a.read(root/label/'config.status.json')['contact_currents_A_per_um']['drain']
        zero=curr('zero');baseline=s.a.read(Path(c['baseline'])/'config.status.json')['contact_currents_A_per_um']['drain'];drift=abs(zero-baseline)
        slopes={amp:float(next(r['derivative_A_per_um'] for r in ports if r['key']==c['key'] and r['port']=='drain' and r['amplitude']==amp)) for amp in ('full','half')}
        linearity=abs(slopes['full']/slopes['half']-1)
        for amp,alpha in (('full',.001),('half',.0005)):
            plus,minus=curr('plus_'+amp),curr('minus_'+amp);odd=(plus-minus)/2;even=abs((plus-zero)+(minus-zero))/2
            error=abs(slopes[amp]/direct-1);snr=abs(odd)/max(drift,1e-300)
            rows.append(dict(key=c['key'],amplitude=amp,tangent_Id_A_per_um=direct,adjoint_Id_A_per_um=dual,duality_relative=dual_error,DC_Id_derivative_A_per_um=slopes[amp],prediction_relative=error,two_amplitude_relative=linearity,even_over_odd=even/abs(odd),signal_over_drift=snr,qualified=error<=1e-3 and dual_error<=1e-6 and linearity<=1e-3 and even/abs(odd)<=.01 and snr>=100))
    s.a.write_csv(OUT/'ports.csv',rows);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'ports.csv']+[p for c in cases for p in (s.LOCAL/'fixed'/c['key']/'zero').glob('adjoint.*')]);print(rows,flush=True)

if __name__=='__main__':run()
