"""Original 816 targets from a cold equilibrium/drain/gate workflow.
No pointwise reclosure: qualify the actual accepted sweep checkpoints.
"""
import argparse,copy,math,os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import simplemos_hfs_cloud_20260926 as h

class Matrix:
    def __init__(self,root,inputs,runner,srh):
        self.root=root.resolve();self.root.mkdir(parents=True,exist_ok=True);self.inputs=inputs.resolve();self.runner=runner.resolve();self.srh=srh.resolve()
        assert h.read(srh/'summary.json')['passed']
        assert h.read(srh/'field_review.json')['reconstruction_passed']
        for p,sha in h.read(inputs/'hashes.json').items():assert h.sha(inputs/p)==sha,p
        self.contract=h.read(inputs/'contract.json');self.executor=object.__new__(h.Run);self.executor.runner=self.runner
        seal=dict(binary_sha256=h.sha(runner),driver_sha256=h.sha(__file__),inputs_sha256=h.sha(inputs/'hashes.json'),srh_gate_sha256=h.sha(srh/'summary.json'),field_gate_sha256=h.sha(srh/'field_review.json'))
        if (root/'seal.json').exists():assert h.read(root/'seal.json')==seal
        else:h.write(root/'seal.json',seal)
    def config(self,c):
        cfg=h.read(self.inputs/c['device']/'template.json')
        for key in ('mesh_file','node_doping_file','materials_file'):cfg[key]=str(self.inputs/cfg[key])
        cfg['simulation_type']='dc_sweep';return cfg
    def sweep(self,c,label,contact,targets,seed=None):
        dest=self.root/c['case']/label;dest.mkdir(parents=True,exist_ok=True);cfg=self.config(c)
        vd=c['vd'] if label=='gate' else 0.
        for x in cfg['contacts']:x['bias']=vd if x['name']=='drain' else 0.
        span=targets[-1]-targets[0]
        cfg['output_csv']=str(dest/'curve.csv')
        if label=='poisson':
            # poisson_block creates only one correction. With that explicit
            # initial state, poisson_only performs the full nonlinear solve.
            # Avoid stopping at a relative reduction of the large cold residual;
            # the following coupled equilibrium retains its original gates.
            cfg['solver']['method']='poisson_only'
            cfg['solver']['reltol']=0.
        sweep=dict(mode='iv',contact=contact,current_contact='drain',start=targets[0],stop=targets[-1],step=.05 if label=='gate' else max(.005,span*.1),bias_points=targets,warm_start=True,write_vtk=False,write_state_file=str(dest/'accepted.h5'),write_state_every_point_prefix=str(dest/'state'),stop_on_failure=True)
        if seed:sweep['initial_state_file']=str(seed)
        else:sweep['initialization']=dict(mode='poisson_block',write_state_file=str(dest/'poisson_handoff.h5'),diagnostic_csv=str(dest/'initialization.csv'))
        if span:
            sweep.update(initial_step=.025 if label=='gate' else span*.1,min_step=span*1e-5,max_step=.125 if label=='gate' else span,growth_factor=1.5,shrink_factor=.5,max_retries=30,continuation=dict(predictor=dict(mode="linear",fields=["psi","phin","phip"],max_extrapolation_ratio=3.)))
        sweep['diagnostics']=dict(newton_history=dict(enabled=True,attempts_csv_file=str(dest/'attempts.csv'),iterations_csv_file=str(dest/'iterations.csv')))
        cfg['sweep']=sweep;h.write(dest/'config.json',cfg);status=self.executor.execute(dest/'config.json')
        print(c['case'],label,'exit',status['exit_code'],'converged',status.get('converged'),flush=True)
        assert status['exit_code']==0 and status.get('converged'),(c,label,status)
        return dest
    def qualify(self,c,dest):
        geo=h.read(self.inputs/c['device']/'geometry.json');mesh=h.read(self.inputs/c['device']/'mesh.json');reference={round(float(r['gate_voltage_V'])/.05):float(r['drain_total_current_A_per_um']) for r in h.rows(self.inputs/c['reference'])};out=[]
        sweep=h.rows(dest/'curve.csv');assert len(sweep)==51
        for r in sweep:
            vg=float(r['bias_V']);index=round(vg/.05);assert abs(vg-index*.05)<1e-12
            state=dest/('state_bias_'+format(vg,'.6f').replace('.','p')+'.h5')
            assert state.exists(),state
            target=dest/f'audit_{index:03d}';target.mkdir(exist_ok=True)
            cfg=self.config(c);cfg.pop('sweep',None)
            for x in cfg['contacts']:x['bias']=c['vd'] if x['name']=='drain' else vg if x['name']=='gate' else 0.
            cfg.update(state_file=str(state),simulation_type='newton_carrier_term_probe',output_csv=str(target/'all_row.csv'),carrier_term_probe=dict(solved_equation_terms=True))
            cfg['solver']['carrier_row_convergence']['mode']='report';cfg['solver']['stable_merit_comparison']=False
            cfg['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
            edge=copy.deepcopy(cfg);edge.pop('carrier_term_probe');edge.update(simulation_type='sg_edge_flux_probe',output_csv=str(target/'acceptance_edges.csv'))
            fun=copy.deepcopy(edge);fun.pop('output_csv');fun.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(target/'port_residual.csv'),contact_edge_output_csv=str(target/'port_edges.csv'))
            for name,deck in [('all_row',cfg),('acceptance_edges',edge),('functional',fun)]:
                h.write(target/(name+'.json'),deck);self.executor.execute(target/(name+'.json'))
            edges=h.rows(target/'acceptance_edges.csv');ports={}
            for contact in mesh['contacts']:
                nodes=set(contact['node_ids'])
                raw={car:1.602176634e-19*1e-6*math.fsum((int(int(e['node0']) in nodes)-int(int(e['node1']) in nodes))*float(e[car+'_particle_line_flux_per_m_s']) for e in edges) for car in ('electron','hole')}
                ports[contact['name']]=-raw['electron']+raw['hole']
            current=float(r['current_total_A_per_um']);point=dict(exit_code=0,converged=r['converged']=='1',contact_currents_A_per_um=ports,iterations=int(r['iterations']),origin='Unmodified accepted dc_sweep checkpoint; convergence from sweep CSV, contacts independently reconstructed from edges',sweep_current_A_per_um=current,state_sha256=h.sha(state))
            h.write(target/'config.status.json',point)
            result=h.qualify(target,geo);diff=abs(ports['drain']/current-1) if current else math.inf
            result.update(**c,index=index,vg=vg,native_Id_A_per_um=reference[index],Id_error_percent=100*(current/reference[index]-1),sweep_port_relative=diff)
            result['qualified']=result['qualified'] and diff<=1e-8
            result['current_accepted']=abs(result['Id_error_percent'])<=2.
            h.write(target/'result.json',result);out.append(result)
            print(c['case'],vg,'qualified',result['qualified'],'Id error %',result['Id_error_percent'],flush=True)
        return out
    def case(self,c):
        try:
            poisson=self.sweep(c,'poisson','gate',[0.])
            eq=self.sweep(c,'equilibrium','gate',[0.],poisson/'accepted.h5')
            drain=self.sweep(c,'drain','drain',[0.,c['vd']],eq/'accepted.h5')
            gate=self.sweep(c,'gate','gate',[i*.05 for i in range(51)],drain/'accepted.h5')
            return self.qualify(c,gate)
        except Exception as e:
            h.write(self.root/c['case']/'failure.json',dict(case=c['case'],error=repr(e)));print('FAILED',c['case'],repr(e),flush=True);return []
    def all(self,pilot):
        cases=[c for c in self.contract['cases'] if not pilot or c['device']=='n23' and c['vd']==1.]
        with ThreadPoolExecutor(max_workers=min(4,len(cases))) as pool:out=[r for group in pool.map(self.case,cases) for r in group]
        h.csvout(self.root/'comparison.csv',[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in out])
        result=dict(requested_points=51*len(cases),completed_points=len(out),qualified_points=sum(r['qualified'] for r in out),current_accepted_points=sum(r['current_accepted'] for r in out),complete=len(out)==51*len(cases) and all(r['qualified'] and r['current_accepted'] for r in out),pilot=pilot)
        h.write(self.root/'summary.json',result);print(result,flush=True)
        assert result['complete'],'Original cold-start matrix incomplete or failed; inspect retained evidence'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--runner',type=Path,required=True);p.add_argument('--srh',type=Path,required=True);p.add_argument('--pilot',action='store_true');a=p.parse_args();Matrix(a.root,a.inputs,a.runner,a.srh).all(a.pilot)
