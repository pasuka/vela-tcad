"""Executable integration checks for the isolated finite-source input contract."""
from pathlib import Path
import copy
import validate_simplemos_srh_finite_20260908 as t

def main():
    t.a.verify(t.OUT/'build_evidence.json');case=t.a.read(t.OUT/'contract.json')['cases'][0]
    geo,mask=t.v.m.previous.prior.support(case);N=geo.count;node=case['mapped_nodes']['792']['node'];contact=int(geo.contact_nodes[0]);root=t.LOCAL/'input_tests'
    valid=f'{N} 1\n{node} 1.1\n';tests=[
        ('wrong_size',f'{N-1} 1\n{node} 1.1\n',{},'Invalid finite SRH node manifest'),
        ('missing_record',f'{N} 1\n',{},'Invalid finite SRH volume record'),
        ('contact_node',f'{N} 1\n{contact} 1.1\n',{},'Invalid finite SRH volume record'),
        ('duplicate',f'{N} 2\n{node} 1.1\n{node} 1.2\n',{},'Invalid finite SRH volume record'),
        ('zero_ratio',f'{N} 1\n{node} 0\n',{},'Invalid finite SRH volume record'),
        ('extra',valid+'unexpected\n',{},'Extra finite SRH data'),
        ('out_of_range',valid,{'VELA_CANDIDATE_SRH_ALPHA':'1.01'},'SRH amplitude outside frozen interpolation range'),
        ('missing_alpha',valid,{'VELA_CANDIDATE_SRH_ALPHA':None},'Missing finite SRH amplitude'),
        ('conflict',valid,{'VELA_CANDIDATE_SRH_NODE':str(node)},'Conflicting SRH diagnostic source selectors')]
    records=[];files=[Path(__file__).resolve(),t.OUT/'build_evidence.json']
    for name,data,changes,expected in tests:
        dest=root/name;dest.mkdir(parents=True,exist_ok=False);path=dest/'volumes.txt';path.write_text(data)
        cfg=t.a.read(t.LOCAL/case['key']/'fixed/zero/functional.json');cfg['output_csv']=str(dest/'unused.csv');cfg['residual_output_csv']=str(dest/'residual.csv');t.a.write(dest/'config.json',cfg)
        env=t.c.env();env.update(VELA_CANDIDATE_SRH_VOLUMES=str(path),VELA_CANDIDATE_SRH_ALPHA='1')
        for k,value in changes.items():
            if value is None:env.pop(k,None)
            else:env[k]=value
        status=t.v.execute(dest/'config.json',t.RUNNER,env);passed=status['exit_code']==1 and expected in status.get('stderr','')
        records.append(dict(name=name,expected=expected,exit_code=status['exit_code'],qualified=passed));files+=list(dest.glob('*'))
    t.a.write_csv(t.OUT/'input_tests.csv',records);t.d.matrix.freeze(t.OUT/'input_tests_evidence.json',files+[t.OUT/'input_tests.csv'])
    assert all(r['qualified'] for r in records),records
    print('Finite SRH invalid input rejection:',len(records),'passed',flush=True)

if __name__=='__main__':main()
