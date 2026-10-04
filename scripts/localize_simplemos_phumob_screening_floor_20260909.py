"""Localize native PhuMob discrepancies to G(P)'s floor, without production edits.

An effective floor inferred from one cell is a diagnostic fit, not an accepted
physical parameter. Held-out cells test whether the discrepancy has this form.
"""
import math
from pathlib import Path
from collections import defaultdict
from decimal import localcontext
import analyze_simplemos_phumob_native_cells_20260909 as c
import audit_simplemos_phumob_cross_derivatives_20260908 as h

a,d=c.a,c.d;OUT=c.OUT/'screening_floor'


def mobility(state,car,G):
    nd,na,n,p=(float(state[k]) for k in ('nd','na','n','p'))
    nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
    mx,mn,nref,alpha,m,other_m=(1417,52.2,9.68e16,.68,1,1.258) if car=='e' else (470.5,44.9,2.23e17,.719,1.258,1)
    other=p if car=='e' else n;nsc=nd+na+other
    P=1/(2.459*nsc**(2/3)/3.97e13+3.828*(n+p)/(1.36e20*m))
    z=P**.6478;F=(.7643*z+2.2999+6.5502*m/other_m)/(z+2.3670-.8552*m/other_m)
    eff=(nd+G*na if car=='e' else na+G*nd)+other/F
    A=mx*mx/(mx-mn)*nsc*(nref/nsc)**alpha+mx*mn/(mx-mn)*(n+p)
    return mx*A/(A+mx*eff)


def main():
    a.verify(c.OUT/'evidence.json')
    a.write(OUT/'contract.json',dict(scope='300 K, same 8 native PhuMob states. No production change and no numerical gate change.',
        classification='Group cell errors by the number of local vertices below the independently computed G minimum.',
        diagnostic_fit='For each carrier infer a single effective G floor from the lexicographically first fully clamped cell. Validate its functional form on every other cell; do not treat fitted agreement as a native parameter or algorithm proof.',
        training_selection='First (key, numeric cell) with all three P below the mathematical minimum and three positive vertex measures.',
        heldout_relative_gate=1e-7))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),c.OUT/'evidence.json',c.OUT/'mapping.json',c.LOCAL/'scalar_input.csv',OUT/'contract.json',Path(h.__file__)])
    inputs={x['id']:x for x in a.rows(c.LOCAL/'scalar_input.csv')};cpp={x['id']:x for x in a.rows(c.LOCAL/'scalar_output.csv')}
    mapping={(x['key'],str(x['cell'])):x for x in a.read(c.OUT/'mapping.json')['records'] if x['mode']=='vertex_local'}
    rows=[r for r in a.rows(c.OUT/'cell_comparison.csv') if r['mode']=='vertex_local']
    with localcontext() as ctx:
        ctx.prec=100
        minima={car:tuple(map(float,h.gminimum(i))) for i,car in enumerate(('e','h'))}
    groups=defaultdict(list);validation=[];fits=[]
    for car in ('e','h'):
        data=sorted([r for r in rows if r['carrier']==car],key=lambda r:(r['key'],int(r['cell'])))
        def count(row):return sum(float(cpp[i]['P_'+car])<minima[car][0] for i in mapping[row['key'],row['cell']]['ids'])
        training=next(r for r in data if count(r)==3 and min(mapping[r['key'],r['cell']]['weights'])>0)
        def prediction(row,floor):
            mp=mapping[row['key'],row['cell']]
            return math.fsum(w*mobility(inputs[i],car,floor if float(cpp[i]['P_'+car])<minima[car][0] else float(cpp[i]['G_'+car])) for i,w in zip(mp['ids'],mp['weights']))
        low,high=0.,1.;target=float(training['native_mu'])
        assert prediction(training,low)>target>prediction(training,high)
        for _ in range(80):
            mid=(low+high)/2
            if prediction(training,mid)>target:low=mid
            else:high=mid
        inferred=(low+high)/2
        fits.append(dict(carrier=car,training_key=training['key'],training_cell=training['cell'],mathematical_P_min=minima[car][0],mathematical_G_min=minima[car][1],inferred_effective_G_floor=inferred,G_difference=inferred-minima[car][1],production_parameter_qualified=False))
        for r in data:
            clamped=count(r);groups[car,clamped].append(abs(float(r['signed_relative'])))
            baseline=prediction(r,minima[car][1]);trial=prediction(r,inferred);native=float(r['native_mu'])
            validation.append(dict(key=r['key'],carrier=car,cell=r['cell'],clamped_vertices=clamped,is_training=r is training,
                independent_formula_relative=baseline/float(r['formula_mu'])-1,original_relative=float(r['signed_relative']),diagnostic_floor_relative=trial/native-1))
    summary=[]
    for (car,num),values in sorted(groups.items()):
        heldout=[r for r in validation if r['carrier']==car and r['clamped_vertices']==num and not r['is_training']]
        summary.append(dict(carrier=car,clamped_vertices=num,cells=len(values),original_max_relative=max(values),heldout_cells=len(heldout),diagnostic_max_relative=max(abs(r['diagnostic_floor_relative']) for r in heldout)))
    a.write_csv(OUT/'floor_inference.csv',fits);a.write_csv(OUT/'cells.csv',validation);a.write_csv(OUT/'groups.csv',summary)
    result=dict(original_cpp_independent_formula_max_relative=max(abs(r['independent_formula_relative']) for r in validation),heldout_diagnostic_max_relative=max(abs(r['diagnostic_floor_relative']) for r in validation if not r['is_training']),
        inferred_constants_are_diagnostic_only=True,native_floor_algorithm_identified=False,production_changed=False)
    a.write(OUT/'summary.json',result)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json']+[OUT/n for n in ('floor_inference.csv','cells.csv','groups.csv','summary.json')])
    print(fits,result,flush=True)


if __name__=='__main__':main()
