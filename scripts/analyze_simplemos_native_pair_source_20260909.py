"""Compare direct Sentaurus DC responses with the frozen Vela pair source."""
import math
import re
from pathlib import Path
import native_simplemos_phumob_pair_source_20260909 as n

a,d=n.a,n.d
OUT=n.OUT/'comparison'
Q_NATIVE=1.602192e-19
Q_VELA=1.602176634e-19

def export_quantum(value):
    """Absolute last-place unit in the native %.14E current export."""
    return 10.**(math.floor(math.log10(abs(value)))-14) if value else 0.

def run():
    a.verify(n.OUT/'response_evidence.json')
    vela=n.s.OUT/'scaled_source'
    a.verify(vela/'port_calibration/evidence.json')
    native=a.rows(n.OUT/'response_points.csv');local=a.rows(vela/'port_calibration/ports.csv');scope=a.rows(n.s.OUT/'scope.csv')
    cases=sorted({r['key'] for r in native})
    a.write(OUT/'contract.json',dict(scope='Four matched integrated pair-source directions, direct native DC current, two signed amplitudes, zero drift and comparison to Vela finite differences. No fit and no finite volume replacement.',
        gates=a.read(n.OUT/'contract.json')['gates'],
        raw_comparison_gate=1e-3,
        charge_normalized_comparison='Supplementary dId/q comparison using the rounded manual q and Vela SI q; the manual printed value is not proof of the native runtime constant. Raw comparison remains independently gated.',
        export_precision='Native current is printed with 15 significant digits. Bound odd-signal rounding by one half of the sum of the two half-last-place errors; use export resolution as a floor on zero-drift noise. Require rounding/signal <= prediction gate.',
        q_native_C=Q_NATIVE,q_vela_C=Q_VELA,
        field_limit='Previous Vela whole-field source calibration remains 6/8; native terminal qualification cannot override it.'))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),OUT/'contract.json',n.OUT/'response_evidence.json',vela/'port_calibration/evidence.json',n.s.OUT/'scope.csv'])
    rows=[];terminals=[];sites=[]
    for key in cases:
        points={r['label']:r for r in native if r['key']==key};assert len(points)==5
        for label,p in points.items():
            raw=n.LOCAL/'response_raw/bundle'/key/label
            exported=(raw/'native_des.plt').read_text().split('Data {',1)[1].split('}',1)[0].split()
            assert exported and all(re.fullmatch(r'-?\d\.\d{14}E[+-]\d+',x) for x in exported), 'Unexpected current export precision'
            seen=[x.split(',') for x in (raw/'pair_source_seen.csv').read_text().splitlines()]
            expected=[r for r in scope if r['key']==key]
            for idx in range(8):
                found=[x for x in seen if int(x[0])==idx];assert found
                r=expected[idx]
                assert all(abs(float(x[1])-float(r['x_um']))<1e-10 and abs(float(x[2])-float(r['y_um']))<1e-10 and float(x[3])==float(r['native_rate_cm3_s']) and float(x[4])==float(p['alpha']) and int(x[5])==2 for x in found)
                integrated=float(found[0][3])*float(r['native_volume_m2'])*1e6
                wanted=float(r['integrated_pair_source_per_m_s'])
                error=abs(integrated/wanted-1) if wanted else abs(integrated)
                assert error<=1e-14
                sites.append(dict(key=key,label=label,site=idx,node=r['node'],source_integral_relative=error,qualified=True))
        derivatives={}
        zero=float(points['zero']['Id_A_per_um'])
        ref=next(r for r in a.rows(n.s.prior.NOUT/'native_points.csv') if r['model']=='phumob' and r['case']==key and int(r['index'])==0)
        drift=abs(zero-float(ref['Id_A_per_um']))
        for amp,alpha in (('full',.001),('half',.0005)):
            for port in ('drain','source','gate','substrate'):
                plus=float(points['plus_'+amp][port+'_A_per_um']);minus=float(points['minus_'+amp][port+'_A_per_um'])
                derivative=(plus-minus)/(2*alpha);terminals.append(dict(key=key,amplitude=amp,port=port,derivative_A_per_um=derivative))
                if port=='drain':derivatives[amp]=derivative
        linearity=abs(derivatives['full']/derivatives['half']-1)
        for amp,alpha in (('full',.001),('half',.0005)):
            plus=float(points['plus_'+amp]['Id_A_per_um']);minus=float(points['minus_'+amp]['Id_A_per_um']);odd=(plus-minus)/2
            even=abs(math.fsum((plus-zero,minus-zero)))/2
            rounding=(export_quantum(plus)+export_quantum(minus))/4
            noise=max(drift,export_quantum(zero),rounding)
            snr=abs(odd)/noise
            rounding_relative=rounding/abs(odd)
            v=next(r for r in local if r['key']==key and r['amplitude']==amp)
            vd=float(v['DC_Id_derivative_A_per_um']);nd=derivatives[amp]
            raw_error=vd/nd-1;charge_error=(vd/Q_VELA)/(nd/Q_NATIVE)-1
            qualified=all(p['qualified']=='True' for p in points.values()) and linearity<=1e-3 and even/abs(odd)<=.01 and snr>=100 and rounding_relative<=1e-3
            rows.append(dict(key=key,amplitude=amp,native_dId_A_per_um=nd,Vela_dId_A_per_um=vd,Vela_over_native_response_minus_one=raw_error,charge_normalized_relative=charge_error,native_two_amplitude_relative=linearity,native_even_over_odd=even/abs(odd),native_signal_over_zero_drift=snr,native_export_rounding_over_signal=rounding_relative,native_response_qualified=qualified,cross_solver_qualified=qualified and v['qualified']=='True' and abs(raw_error)<=1e-3))
    a.write_csv(OUT/'responses.csv',rows);a.write_csv(OUT/'terminal_derivatives.csv',terminals);a.write_csv(OUT/'source_identity.csv',sites)
    summary=dict(native_DC=len(native),qualified_native_DC=sum(r['qualified']=='True' for r in native),source_site_checks=len(sites),qualified_source_sites=len(sites),native_responses=len(rows),qualified_native_responses=sum(r['native_response_qualified'] for r in rows),qualified_cross_solver_responses=sum(r['cross_solver_qualified'] for r in rows),maximum_raw_response_relative=max(abs(r['Vela_over_native_response_minus_one']) for r in rows),maximum_charge_normalized_relative=max(abs(r['charge_normalized_relative']) for r in rows),previous_Vela_field_response_qualified='6/8',finite_volume_replacement=False)
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'responses.csv',OUT/'terminal_derivatives.csv',OUT/'source_identity.csv',OUT/'summary.json'])
    print(summary,flush=True)

if __name__=='__main__':run()
