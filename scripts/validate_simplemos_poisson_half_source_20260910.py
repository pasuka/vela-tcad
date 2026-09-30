"""Complete the frozen second source amplitude for the packed Poisson candidate."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_poisson_precision_20260910 as c
s=c.s;LOCAL=c.LOCAL/'half_amplitude';OUT=c.OUT/'half_amplitude'

def run():
    s.a.verify(c.OUT/'build_evidence.json');s.a.verify(c.OUT/'kernel_check/evidence.json')
    jobs=[];files=[Path(__file__).resolve(),c.OUT/'build_evidence.json',c.OUT/'kernel_check/evidence.json',c.m.OUT/'dc_evidence.json']
    for case in s.a.read(c.m.OUT/'contract.json')['jobs']:
        if case['kind']!='source' or case['label'] not in ('plus_half','minus_half'):continue
        base=Path(case['dest']);dest=LOCAL/case['key']/case['label'];cfg=s.a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv')
        s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);jobs.append(dict(case,base=str(base),dest=str(dest),mode='packed'));files+=list(dest.glob('*.json'))+[Path(cfg['state_file']),Path(case['source_file'])]
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Same packed Poisson candidate, source values, original settled bases, gates and 200-iteration limit; add the previously frozen +/-0.0005 amplitude. Full-amplitude results are required separately; this contract does not presume they pass.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    def one(j):
        dest=Path(j['dest']);env=s.env(j,j['alpha']);env.update(VELA_STABLE_MERIT='1',VELA_POISSON_PRECISION='packed',VELA_STABLE_MERIT_TRACE=str(dest/'merit.csv'))
        s.V.execute(dest/'config.json',c.RUNNER,env)
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),c.RUNNER,env)
        result=dict(**j,**s.prior.q.run.w.old.prior.old.qualify(j,dest));print(result,flush=True);return result
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,jobs))
    s.a.write_csv(OUT/'dc.csv',rows);s.d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'freeze.json',OUT/'dc.csv']+[p for j in jobs for p in Path(j['dest']).rglob('*') if p.is_file()])

if __name__=='__main__':run()
