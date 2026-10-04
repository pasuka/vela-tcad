"""Check whether a Poisson-only split extension preserves full-DD partition invariance."""
import math,subprocess
from decimal import Decimal as D,localcontext
from pathlib import Path
import validate_simplemos_split_psi_v2_20260911 as v

def main():
    v.verify(v.OUT/'replay_evidence.json');dest=v.LOCAL/'partition';dest.mkdir(parents=True,exist_ok=False)
    cpp=v.CPP.replace('if(mode=="split") {nn[i].xpsi=pair.hi;low[i]=pair.lo;}',
       'if(mode=="split") {nn[i].xpsi=pair.hi;low[i]=pair.lo; if(t.contains("low")){nn[i].xpsi=t["candidate"][i];low[i]=t["low"][i];}}')
    (dest/'kernel.cpp').write_text(cpp);cmd=['D:/msys64/ucrt64/bin/c++.exe','-std=c++20','-O2','-I'+str(v.ROOT/'include'),str(dest/'kernel.cpp'),'-o',str(dest/'kernel.exe')]
    v.write(dest/'compile.json',cmd);p=subprocess.run(cmd,env=v.env(),capture_output=True,text=True);(dest/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    files=[Path(__file__).resolve(),v.OUT/'replay_evidence.json'];records=[]
    for c in v.rows(v.OUT/'identity.csv'):
        base=Path(c['new_dest']);d=dest/c['device']/c['vd']/c['index'];cfg=v.read(base/'config.json');seed=Path(cfg['state_file']);state=v.rows(seed)
        data=v.read(base/'inputs.json');trial=data['trials'][0];n=len(state);bc={int(i) for i,_ in data['bcpsi']};low=[0.]*n;candidate=list(trial['candidate'])
        with localcontext() as ctx:
            ctx.prec=100
            for i,r in enumerate(state):
                original=float(r['packed_psi']);assert original==trial['x'][i]
                if i not in bc:
                    hi=math.nextafter(original,math.inf);low[i]=original-hi
                    assert D.from_float(hi)+D.from_float(low[i])==D.from_float(original)
                    r['packed_psi']=format(hi,'.17g');r['psi']=format(hi*float(r['packed_potential_scale_V']),'.17g');candidate[i]=hi
        v.csvout(d/'state.csv',state);cfg.update(state_file=str(d/'state.csv'),residual_output_csv=str(d/'residual.csv'),contact_edge_output_csv=str(d/'contact.csv'));v.write(d/'config.json',cfg)
        assert v.execute(d/'config.json',v.ROOT/'build-release/vela_example_runner.exe')==0
        alt=dict(trial,label='partition',candidate=candidate,low=low);data['trials']=[trial,alt];v.write(d/'inputs.json',data)
        p=subprocess.run([str(dest/'kernel.exe'),str(d/'inputs.json'),str(d/'kernel.csv')],env=v.env(),capture_output=True,text=True);assert p.returncode==0,p.stderr
        kr=v.rows(d/'kernel.csv');original=[r['R'] for r in kr if r['trial']=='base' and r['mode']=='split'];alt=[r['R'] for r in kr if r['trial']=='partition' and r['mode']=='split'];assert original==alt
        baseR=v.rows(base/'residual.csv');otherR=v.rows(d/'residual.csv');delta={}
        for name in ('phin_residual','phip_residual'):
            diff=[float(y[name])-float(x[name]) for x,y in zip(baseR,otherR)]
            delta[name+'_changed']=sum(x!=0 for x in diff);delta[name+'_max_difference']=max(map(abs,diff))
        rec=dict(device=c['device'],vd=c['vd'],index=c['index'],same_exact_potential=True,poisson_identical=True,**delta)
        rec['full_partition_invariant']=delta['phin_residual_changed']==delta['phip_residual_changed']==0
        records.append(rec);print(rec,flush=True)
        files += [seed]+[p for p in d.iterdir() if p.is_file()]
    v.csvout(v.OUT/'partition.csv',records)
    v.freeze(v.OUT/'partition_evidence.json',files+[v.OUT/'partition.csv']+[p for p in dest.iterdir() if p.is_file()])

if __name__=='__main__':main()
