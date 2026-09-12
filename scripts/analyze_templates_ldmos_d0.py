"""Score full D0 electrical curves and all native nodal temperature fields.

Uses unchanged Stage-4 electrical limits and the separately approved D0 thermal
contract. A zero-power point retains the undefined relative heat balance, and
must independently prove exactly zero heat and constant ambient temperature.
"""
import argparse,hashlib,json
from pathlib import Path
from analyze_templates_ldmos_stage4_d5 import curve_error,ratio_error,kcl_audit,verdict,read_curve,read_vela_curve
from analyze_templates_ldmos_thermal import assess
from run_templates_ldmos_electrothermal_curve import state_gate


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def analyze(native,curves,contract):
    manifest=read(native/'manifest.json');fields=manifest['fields'];reference={};candidate={};metrics={};thermal=[];sources={};numerical=[];continuation=[]
    for gate,index in ((4,1),(8,2)):
        key=f'Vg{gate}';directory=curves[gate];ledger=read(directory/'ledger.json')
        if ledger['status']!='complete' or len(ledger['exact_points'])!=31:raise ValueError('Full completed 31-point curves required')
        if not ledger.get('runs'):raise ValueError('Continuation evidence is required')
        for run in ledger['runs']:
            if not run['gate']['pass_gate']:continue
            path=Path(run['directory'])/'output.json';check=state_gate(read(path),run['bias_V'])
            continuation.append(dict(gate=gate,bias_V=run['bias_V'],**check));sources[str(path)]=sha(path)
        refpath=native/'normalized'/f'IdVd_Vg{index}_n4_des_drain_curve.csv'
        reference[key]=read_curve(refpath);candidate[key]=read_vela_curve(directory/'curve.csv')
        if len(reference[key])!=31 or len(candidate[key])!=31:raise ValueError('Exact 31-point current curves required')
        metrics[key]=curve_error(reference[key],candidate[key])
        selected=sorted((f for f in fields if f['gate_V']==gate),key=lambda f:f['point_index'])
        if [f['point_index'] for f in selected]!=list(range(31)):raise ValueError('Missing/duplicate native temperature field')
        for point,field in zip(ledger['exact_points'],selected):
            if abs(point['bias_V']-field['bias_V'])>1e-9:raise ValueError('Native temperature bias mismatch')
            rp=Path(point['result']);tp=Path(field['temperature_file']);result=read(rp);ref=read(tp)
            if sha(tp)!=field['temperature_sha256']:raise ValueError('Native temperature hash mismatch')
            if ref['node_id']!=list(range(len(result['temperature_K']))):raise ValueError('Temperature node mapping mismatch')
            check=state_gate(result,point['bias_V']);numerical.append(dict(gate=gate,bias_V=point['bias_V'],**check))
            score=assess(ref['temperature_K'],result['temperature_K'],result['nodal_area_m2'],result['lattice_source_W_per_m'],result['boundary_heat_W_per_m'],contract)
            score.update(gate=gate,bias_V=point['bias_V'])
            if point['bias_V']==0.:
                score['zero_power_equilibrium_pass']=check['pass_gate'] and result['lattice_source_W_per_m']==0. and result['boundary_heat_W_per_m']==0. and all(t==contract['ambient_temperature_K'] for t in result['temperature_K'])
                score['qualification_pass']=score['zero_power_equilibrium_pass'] and all(v['status']=='pass' for k,v in score['gates'].items() if k!='heat_balance')
            else:score['qualification_pass']=score['status']=='pass'
            thermal.append(score);sources[str(rp)]=sha(rp);sources[str(tp)]=sha(tp)
        for p in (refpath,directory/'curve.csv',directory/'terminal_balance.csv',directory/'ledger.json'):sources[str(p)]=sha(p)
    metrics['gate_ratio_error_percent']=ratio_error(reference,candidate)
    metrics['kcl_audit']={f'Vg{g}':kcl_audit(curves[g]/'terminal_balance.csv') for g in (4,8)}
    metrics['max_normalized_kcl_percent']=max(v['max_normalized_kcl_percent'] for v in metrics['kcl_audit'].values())
    electrical=verdict(metrics,'final');passed=electrical['status']=='pass' and all(t['qualification_pass'] for t in thermal) and all(n['pass_gate'] for n in numerical+continuation)
    return dict(schema='vela.templates_ldmos.d0_summary.v1',status='pass' if passed else 'fail',scope='Experimental four-equation SI adapter; ideal contact hRec equivalence retains its separately bounded scope; no Thermodynamic, Peltier or RecGenHeat',electrical=electrical,metrics=metrics,thermal=thermal,numerical=numerical,continuation=continuation,thermal_contract=contract,sources_sha256=sources)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('native','vg4','vg8','contract','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=analyze(a.native,{4:a.vg4,8:a.vg8},read(a.contract))
    result['sources_sha256'][str(a.contract)]=sha(a.contract);result['sources_sha256'][str(a.native/'manifest.json')]=sha(a.native/'manifest.json')
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(result['status'])
if __name__=='__main__':main()
