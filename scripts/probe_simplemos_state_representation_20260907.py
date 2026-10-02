"""Representation-only follow-up to retained native-state port mismatches."""
import copy
from decimal import Decimal,localcontext
from pathlib import Path
import prepare_simplemos_production_consistency_20260907 as audit

a=audit.a; d=audit.d; p=audit.p; OUT=audit.OUT; LOCAL=audit.LOCAL


def main():
    a.verify(OUT/'probe_evidence.json'); cases=a.read(OUT/'contract.json')['cases']; files=[Path(__file__).resolve(),OUT/'probe_evidence.json']; errors=[]
    jobs=[]
    for c in cases:
        root=LOCAL/c['key']/'native_referenced_joint'; rows=a.rows(Path(c['native_initial']))
        refs=a.rows(Path(c['root'])/'joint/replacement/state.csv'); assert len(rows)==len(refs)
        maximum=0.
        with localcontext() as ctx:
            ctx.prec=70
            for s,r in zip(rows,refs):
                assert s['node_id']==r['node_id']
                for car,field in [('electron','phin'),('hole','phip')]:
                    original=Decimal.from_float(float(s[field])); ref=float(r[car+'_qf_reference_V']); inc=float(original-Decimal.from_float(ref))
                    error=abs(Decimal.from_float(ref)+Decimal.from_float(inc)-original); maximum=max(maximum,float(error))
                    s[car+'_qf_reference_V']=format(ref,'.17g'); s[car+'_qf_increment_V']=format(inc,'.17g')
        assert maximum<=2e-16
        a.write_csv(root/'state.csv',rows); files.append(root/'state.csv'); errors.append(dict(key=c['key'],maximum_binary_physical_qf_change_V=maximum))
        for name in ('functional','edges','terms'):
            old=LOCAL/c['key']/'native_state_joint'/(name+'.json'); cfg=a.read(old); cfg['state_file']=str(root/'state.csv')
            for key in ('output_csv','residual_output_csv','contact_edge_output_csv'):
                if key in cfg:cfg[key]=str(root/Path(cfg[key]).name)
            path=root/(name+'.json'); a.write(path,cfg); files.append(path)
            jobs.append(dict(case=c,config=str(path),strength=1.))
    a.write(OUT/'representation_contract.json',dict(status='frozen_before_representation_control',jobs=jobs,qf_binary_change_limit_V=2e-16,
        explanation='Same native-coherent state stored in the existing contact-basin reference/increment form. Psi and densities unchanged; no solver, coefficient or physical-model change. Original unreferenced port failures retained.'))
    a.write_csv(OUT/'representation_input_changes.csv',errors); files += [OUT/'representation_contract.json',OUT/'representation_input_changes.csv']
    d.matrix.freeze(OUT/'representation_freeze.json',files)
    import subprocess,json,time
    for job in jobs:
        cfg=Path(job['config']); c=job['case']; env=audit.b.env(cfg.parent,Path(c['ratios']),1.,1.)
        for key in list(env):
            if key.startswith(('VELA_DIAGNOSTIC_','VELA_VALIDATE_','VELA_MINORITY_','VELA_SIMPLEMOS_LINEAR_')):env.pop(key)
        start=time.monotonic(); r=subprocess.run([str(audit.b.RUNNER),'--config',str(cfg),'--log','off'],env=env,capture_output=True,text=True)
        cfg.with_suffix('.stdout.txt').write_text(r.stdout); cfg.with_suffix('.stderr.txt').write_text(r.stderr)
        s=json.loads(r.stdout.strip().splitlines()[-1]); s.update(exit_code=r.returncode,elapsed_seconds=time.monotonic()-start); a.write(cfg.with_suffix('.status.json'),s); assert r.returncode==0
    result=[]
    for c in cases:
        old=a.read(LOCAL/c['key']/'native_state_joint/functional.status.json'); new=a.read(LOCAL/c['key']/'native_referenced_joint/functional.status.json')
        before=abs(old['current_A_per_um']/old['contact_current_extractor_A_per_um']-1); after=abs(new['current_A_per_um']/new['contact_current_extractor_A_per_um']-1)
        result.append(dict(key=c['key'],original_relative=before,referenced_relative=after,original_qualified=before<=1e-8,referenced_qualified=after<=1e-8,
            functional_change_A_per_um=new['current_A_per_um']-old['current_A_per_um'],extractor_change_A_per_um=new['contact_current_extractor_A_per_um']-old['contact_current_extractor_A_per_um']))
    a.write_csv(OUT/'representation_port_comparison.csv',result)
    files += [x for c in cases for x in (LOCAL/c['key']/'native_referenced_joint').rglob('*') if x.is_file()]+[OUT/'representation_port_comparison.csv']
    d.matrix.freeze(OUT/'representation_evidence.json',files)
    print(result,flush=True)


if __name__=='__main__':main()
