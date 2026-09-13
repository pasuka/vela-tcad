"""Summarize bounded D0 representative-point and local-field comparisons.

This does not apply full-curve electrical gates to a four-point subset.
"""
import argparse,hashlib,json
from pathlib import Path
from analyze_templates_ldmos_stage4_d5 import read_curve

def read(path):return json.loads(path.read_text(encoding='utf-8'))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('points','local','legacy-local','native','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();points=read(a.points);local=read(a.local);legacy=read(a.legacy_local)
    if points['status']!='pass' or local['status']!='completed' or legacy['status']!='completed':raise ValueError('Incomplete evidence')
    expected={(g,i) for g in (4,8) for i in (0,1,10,30)}
    if {(r['gate'],r['index']) for r in points['points']}!=expected or len(points['points'])!=8:raise ValueError('Expected eight representative points')
    for source in (local,legacy):
        if {(c['gate_V'],c['index']) for c in source['cases']}!=expected or len(source['cases'])!=8:raise ValueError('Incomplete local audit')
    sources=[a.points,a.local,a.legacy_local,Path(__file__)];electrical=[];per_gate={}
    for gate,tag in ((4,1),(8,2)):
        path=a.native/'normalized'/f'IdVd_Vg{tag}_n4_des_drain_curve.csv';sources.append(path);reference=read_curve(path)
        if len(reference)!=31:raise ValueError('Expected native 31-point curve after documented duplicate-zero handling')
        rows=[r for r in points['points'] if r['gate']==gate]
        for row in rows:
            matches=[current for bias,current in reference if abs(bias-row['bias_V'])<1e-9]
            if len(matches)!=1:raise ValueError('Ambiguous native voltage match')
            native=matches[0];candidate=row['current_A_per_m']*1e-6
            electrical.append(dict(gate_V=gate,index=row['index'],bias_V=row['bias_V'],native_A_per_um=native,vela_A_per_um=candidate,
                                   relative_error_percent=100*abs(candidate/native-1) if row['bias_V'] else None))
        per_gate[gate]=dict(reclosure_updates=sum(r['updates'] for r in rows),wall_seconds=sum(r['wall_seconds'] for r in rows),
                           max_selected_current_error_percent=max(r['relative_error_percent'] for r in electrical if r['gate_V']==gate and r['relative_error_percent'] is not None),
                           max_peak_temperature_error_K=max(r['thermal']['gates']['peak_temperature_rise']['error_K'] for r in rows),
                           max_temperature_field_rms_error_K=max(r['thermal']['gates']['temperature_rise_field']['rms_error_K'] for r in rows))
    field=[]
    for row in local['cases']:
        old=next(c for c in legacy['cases'] if c['name']==row['name'])
        field.append(dict(name=row['name'],gate_V=row['gate_V'],bias_V=row['bias_V'],
                          auger_max_abs_before=old['fixed_native']['auger_m3_per_s']['max_abs'],
                          auger_max_abs_after=row['fixed_native']['auger_m3_per_s']['max_abs'],
                          fixed_native=row['fixed_native'],self_consistent=row['reclosed_A'],native_state_max_difference=row['native_state_max_difference']))
    result=dict(status='completed',scope=__doc__,representative_numeric_and_thermal_gates='pass',full_curve_qualification=False,
                per_gate=per_gate,electrical_points=electrical,fields=field,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(per_gate,indent=2))

if __name__=='__main__':main()
