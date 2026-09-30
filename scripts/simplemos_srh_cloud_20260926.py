"""SRH-only signed-Si quadrature controls on the frozen merged HFS baseline.
Large states stay on the compute host. No historical baseline is modified.
"""
import argparse,copy,json,math,os,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal,localcontext
import numpy as np
import simplemos_hfs_cloud_20260926 as h


def signed_si(mesh):
    nodes={r['id']:(r['x']*1e-6,r['y']*1e-6) for r in mesh['nodes']}
    silicon={r['id'] for r in mesh['regions'] if r['material']=='Si'}
    areas=np.zeros(len(nodes))
    for t in mesh['triangles']:
        if t['region_id'] not in silicon:continue
        ids=t['node_ids'];p=[nodes[i] for i in ids];cot=[];l=[]
        for k in range(3):
            a,b,c=p[k],p[(k+1)%3],p[(k+2)%3]
            u=(b[0]-a[0],b[1]-a[1]);v=(c[0]-a[0],c[1]-a[1])
            cot.append((u[0]*v[0]+u[1]*v[1])/abs(u[0]*v[1]-u[1]*v[0]))
            l.append(u[0]*u[0]+u[1]*u[1])
        for k,i in enumerate(ids):areas[i]+=(l[k]*cot[(k+2)%3]+l[(k+2)%3]*cot[(k+1)%3])/8
    return areas


def state_vector(result,geo):
    rr=h.ordered(Path(result['dest'])/'state.csv',geo['count']);out=[]
    with localcontext() as ctx:
        ctx.prec=100
        for name,ref in [('psi',None),('electron_qf_increment','electron_qf_reference_V'),('hole_qf_increment','hole_qf_reference_V')]:
            for i in geo['free_si']:
                r=rr[i];x=(Decimal.from_float(float(r['packed_'+name]))+Decimal.from_float(float(r['packed_'+name+'_low'])))*Decimal.from_float(float(r['packed_potential_scale_V']))
                if ref:x+=Decimal.from_float(float(r[ref]))
                out.append(x)
    return out


class Experiment:
    def __init__(self,root,baseline,runner):
        self.root=root.resolve();self.root.mkdir(parents=True,exist_ok=True)
        self.base=baseline.resolve();self.runner=runner.resolve();self.runs={}
        self.cases=h.read(self.base/'inputs/contract.json')['cases']
        self.old=h.read(self.base/'summary/control_selected.json')
        assert h.read(self.base/'summary/progress.json')['numerical_baseline_complete']
        contract=dict(scope='SRH-only signed Si area, 300K merged HFS; native physical models and all numerical gates unchanged',fractions=[0.,.001,-.001,.0005,-.0005,1.],response_gates=dict(two_amplitude_relative=1e-3,even_over_odd=.01),fixed_source_relative=1e-10,jvp_entry_relative=1e-4,Id_acceptance_percent=2.,legacy_Id_threshold_percent=1.,baseline_seal=h.sha(self.base/'run_seal.json'),binary_sha256=h.sha(runner),driver_sha256=h.sha(__file__),note='2% follows user authorization; no relaxation of residual, KCL, dual or reconstruction gates. Current scope only; original matrix requires separate qualification.')
        if (self.root/'contract.json').exists():assert h.read(self.root/'contract.json')==contract
        else:h.write(self.root/'contract.json',contract)
    def run(self,f):
        if f in self.runs:return self.runs[f]
        name={0.:'zero',1.:'finite',.001:'plus_full',-.001:'minus_full',.0005:'plus_half',-.0005:'minus_half'}[f]
        root=self.root/name;(root/'inputs').mkdir(parents=True,exist_ok=True)
        h.write(root/'inputs/contract.json',h.read(self.base/'inputs/contract.json'))
        for c in self.cases:
            src=self.base/'inputs'/c['case'];dst=root/'inputs'/c['case'];dst.mkdir(exist_ok=True)
            for p in src.iterdir():
                if p.name=='template.json':continue
                if not (dst/p.name).exists():(dst/p.name).symlink_to(p,target_is_directory=p.is_dir())
            cfg=h.read(src/'template.json');cfg['solver']['region_resolved_interface_assembly']['srh_signed_transport_volume_fraction']=f
            h.write(dst/'template.json',cfg)
        if not (root/'seeds').exists():(root/'seeds').symlink_to(self.base/'seeds',target_is_directory=True)
        h.write(root/'run_seal.json',dict(contract_sha256=h.sha(self.root/'contract.json'),fraction=f,templates={c['case']:h.sha(root/'inputs'/c['case']/'template.json') for c in self.cases}))
        self.runs[f]=h.Run(root,self.runner,4);return self.runs[f]
    def seed(self,c,index):
        return next(r for r in self.old if r['case']==c['case'] and r['index']==index and r['arm']=='enormal')
    def execute_targets(self,f,indices,dual=False):
        run=self.run(f)
        jobs=[(c,i,arm) for c in self.cases for i in indices for arm in (('baseline','native') if dual else ('baseline',))]
        def one(j):
            c,i,arm=j;seed=Path(self.seed(c,i)['dest'])/'state.h5' if arm=='baseline' else self.base/'seeds'/c['case']/f'native/vg_{i:03d}/state.h5'
            return run.target(c,i,arm,seed,'controls')
        with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,jobs))
        h.write(run.summary/'selected.json',results)
        assert all(r['qualified'] for r in results),f'{f}: numerical qualification failed'
        return results
    def preflight(self):
        zero=self.execute_targets(0.,[0,10,40,50]);checks=[]
        for r in zero:
            old=self.seed(r,r['index']);d=h.dual(r,old,self.run(0.).geo[r['case']]);d.update(case=r['case'],index=r['index'])
            checks.append(d)
        h.csvout(self.root/'zero_identity.csv',checks);assert all(r['qualified'] for r in checks)
        def fixed(c):
            cfg=h.read(Path(self.seed(c,0)['dest'])/'all_row.json');cfg['state_file']=str(Path(self.seed(c,0)['dest'])/'state.h5')
            terms={};records=[]
            for f in (0.,1.):
                run=self.run(f);dest=run.root/'fixed'/c['case'];dest.mkdir(parents=True,exist_ok=True)
                deck=copy.deepcopy(cfg);deck['solver']['region_resolved_interface_assembly']['srh_signed_transport_volume_fraction']=f;deck['output_csv']=str(dest/'terms.csv')
                h.write(dest/'terms.json',deck);assert run.execute(dest/'terms.json')['exit_code']==0
                terms[f]=h.ordered(dest/'terms.csv',run.geo[c['case']]['count'])
            geo=self.run(0.).geo[c['case']];areas=signed_si(h.read(self.base/'inputs'/c['case']/'mesh_file.json'));err=0.
            for i in geo['free_si']:
                ratio=areas[i]/geo['all_cell'][i]
                for car in ('electron','hole'):
                    a,b=terms[0.][i],terms[1.][i]
                    assert a[car+'_flux_abs_sum']==b[car+'_flux_abs_sum']
                    expected=float(a[car+'_recombination'])*ratio;actual=float(b[car+'_recombination'])
                    if expected:err=max(err,abs(actual/expected-1))
            run=self.run(1.);dest=run.root/'fixed'/c['case'];deck=copy.deepcopy(cfg)
            deck['solver']['region_resolved_interface_assembly']['srh_signed_transport_volume_fraction']=1.
            deck.pop('carrier_term_probe',None);deck.update(simulation_type='newton_jvp_probe',output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'jvp_rows.csv'),directions=[dict(name=f'{mode}_{step}',mode=mode,amplitude_V=step,node_ids=[792,795,1056],exclude_contacts=False) for mode in ('psi','phin','phip') for step in (1e-12,1e-20)],sample_rows=[dict(block=b,node_id=i) for b in ('psi','phin','phip') for i in range(geo['count'])])
            h.write(dest/'jvp.json',deck);assert run.execute(dest/'jvp.json')['exit_code']==0
            bad=[];maximum=0.
            for r in h.rows(dest/'jvp_rows.csv'):
                a,b=float(r['analytic_derivative']),float(r['finite_difference_derivative']);scale=max(abs(a),abs(b));rel=abs(a-b)/scale if scale else 0.
                maximum=max(maximum,rel)
                if not math.isfinite(rel) or rel>1e-4:bad.append(dict(**r,true_relative=rel))
            h.write(dest/'jvp_failures.json',bad)
            result=dict(**c,source_relative=err,jvp_max_relative=maximum,jvp_failures=len(bad),qualified=err<=1e-10 and not bad)
            print('preflight',result,flush=True);return result
        # Initialize runs serially before threaded file reads.
        self.run(1.)
        with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fixed,self.cases))
        h.csvout(self.root/'preflight.csv',records);assert all(r['qualified'] for r in records),'Preflight failed'
    def calibration(self):
        results={0.:h.read(self.run(0.).summary/'selected.json')}
        for f in (.001,-.001,.0005,-.0005):results[f]=self.execute_targets(f,[0])
        checks=[]
        for c in self.cases:
            group={f:next(r for r in rr if r['case']==c['case'] and r['index']==0) for f,rr in results.items()}
            geo=self.run(0.).geo[c['case']]
            with localcontext() as ctx:
                ctx.prec=100
                v={f:state_vector(r,geo) for f,r in group.items()};base=v[0.]
                derivative={a:np.array([float((x-y)/Decimal(str(2*a))) for x,y in zip(v[a],v[-a])]) for a in (.001,.0005)}
                full,half=derivative[.001],derivative[.0005]
                norm=np.linalg.norm(half);relative=np.linalg.norm(full-half)/norm if norm else math.inf
                even=np.array([float((x+y-2*z)/2) for x,y,z in zip(v[.001],v[-.001],base)])
                parity=np.linalg.norm(even)/(.001*np.linalg.norm(full)) if norm else math.inf
            idd={a:(group[a]['current_A_per_um']-group[-a]['current_A_per_um'])/(2*a) for a in (.001,.0005)}
            idrel=abs(idd[.001]/idd[.0005]-1) if idd[.0005] else math.inf
            row=dict(**c,state_two_amplitude_relative=relative,even_over_odd=parity,Id_two_amplitude_relative=idrel,Id_derivative_A_per_um=idd[.0005],qualified=relative<=1e-3 and parity<=.01 and idrel<=1e-3)
            checks.append(row);print('response',row,flush=True)
        h.csvout(self.root/'response.csv',checks);assert all(r['qualified'] for r in checks),'Response calibration failed'
    def finite(self):
        results=self.execute_targets(1.,[0,10,40,50],True);checks=[]
        native={(r['case'],int(r['index'])):float(r['Id_A_per_um']) for r in h.rows(self.base/'inputs/native_points.csv')}
        for c in self.cases:
            for i in (0,10,40,50):
                pair=[r for r in results if r['case']==c['case'] and r['index']==i];a=next(r for r in pair if r['arm']=='baseline');old=self.seed(c,i)
                d=h.dual(*pair,self.run(1.).geo[c['case']]);d.update(**c,index=i,vg=i/50,Id_A_per_um=a['current_A_per_um'],Id_error_percent=100*(a['current_A_per_um']/native[c['case'],i]-1),baseline_Id_error_percent=100*(old['current_A_per_um']/native[c['case'],i]-1))
                d['current_accepted']=abs(d['Id_error_percent'])<=2.;checks.append(d)
        h.csvout(self.root/'finite_comparison.csv',checks)
        passed=all(r['qualified'] and r['current_accepted'] for r in checks)
        h.write(self.root/'summary.json',dict(stage='finite_controls_complete',states=len(results),qualified_states=sum(r['qualified'] for r in results),pairs=len(checks),qualified_pairs=sum(r['qualified'] for r in checks),max_abs_Id_error_percent=max(abs(r['Id_error_percent']) for r in checks),passed=passed,next_stage='Original 816-point and cold-start preparation; requires reviewed source/field results'))
        assert passed,'Finite gate failed'
    def all(self):self.preflight();self.calibration();self.finite()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--runner',type=Path,required=True);a=p.parse_args()
    Experiment(a.root,a.baseline,a.runner).all()
