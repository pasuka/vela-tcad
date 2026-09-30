"""Fail-closed input controls for the isolated additive-source hook."""
from pathlib import Path
import validate_simplemos_phumob_pair_source_scaled_20260909 as checked
s=checked.s
LOCAL=s.LOCAL/'input_guards';OUT=s.OUT/'input_guards'

def run():
    c=s.a.read(s.OUT/'contract.json')['cases'][0];valid=Path(c['source_file']).read_text();lines=valid.splitlines();N=int(lines[0].split()[0]);geo,_=s.fields.support(c)
    tests=[('size',valid.replace(str(N)+' 8',str(N-1)+' 8',1),'0','Invalid pair-source manifest'),
           ('count',valid.replace(str(N)+' 8',str(N)+' 7',1),'0','Invalid pair-source manifest'),
           ('duplicate','\n'.join(lines[:-1]+[lines[1]])+'\n','0','Invalid pair-source record'),
           ('contact','\n'.join([lines[0],f'{geo.contact_nodes[0]} 1']+lines[2:])+'\n','0','Invalid pair-source record'),
           ('truncated','\n'.join(lines[:-1])+'\n','0','Invalid pair-source record'),
           ('extra',valid+'unexpected\n','0','Extra pair-source data'),
           ('alpha_range',valid,'1.1','Invalid pair-source amplitude'),
           ('alpha_trailing',valid,'0junk','Invalid pair-source amplitude'),
           ('alpha_missing',valid,None,'Missing pair-source amplitude')]
    files=[];jobs=[]
    for label,text,alpha,error in tests:
        root=LOCAL/label;root.mkdir(parents=True,exist_ok=False);path=root/'source.txt';path.write_text(text)
        cfg=s.a.read(s.LOCAL/'fixed'/c['key']/'zero/terms.json');cfg['output_csv']=str(root/'terms.csv');s.a.write(root/'config.json',cfg)
        jobs.append(dict(label=label,source=str(path),alpha=alpha,error=error));files += [path,root/'config.json']
    s.a.write(OUT/'contract.json',dict(jobs=jobs,expected='All nine rejected before emitting a terms CSV. No change to accepted output or thresholds.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[Path(__file__).resolve(),OUT/'contract.json',s.OUT/'build_evidence.json']);rows=[]
    for job in jobs:
        root=LOCAL/job['label'];e=s.env(c,0.);e['VELA_PAIR_SOURCE_FILE']=job['source']
        if job['alpha'] is None:e.pop('VELA_PAIR_SOURCE_ALPHA')
        else:e['VELA_PAIR_SOURCE_ALPHA']=job['alpha']
        st=s.V.execute(root/'config.json',s.RUNNER,e);ok=st['exit_code']!=0 and job['error'] in st.get('stderr','') and not (root/'terms.csv').exists()
        rows.append(dict(label=job['label'],exit_code=st['exit_code'],qualified=ok));assert ok,st
    s.a.write_csv(OUT/'checks.csv',rows);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'checks.csv']+[p for p in LOCAL.rglob('*') if p.is_file()]);print('Input guards',len(rows),'passed',flush=True)

if __name__=='__main__':run()
