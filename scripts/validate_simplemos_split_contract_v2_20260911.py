"""Export current geometry and qualify the complete experimental split-state operator."""
import argparse,csv,json,math,os,shlex,subprocess
from pathlib import Path
import validate_simplemos_split_psi_v2_20260911 as prior
ROOT=prior.ROOT;LOCAL=ROOT/'build-release/split_contract_0911/v2';OUT=ROOT/'reference_tcad/simplemos_sentaurus2022/phumob_split_contract_20260911/v2'
read,rows,write,freeze,verify,env,csvout=prior.read,prior.rows,prior.write,prior.freeze,prior.verify,prior.env,prior.csvout
HOOK=r'''
        if(const char* path=std::getenv("VELA_CONTRACT_EXPORT")) {
            static bool written=false;
            if(!written && !bcs.psi.empty()) {
                written=true;
                if(mobilityConfig_.model!="phumob" || mobilityConfig_.edgeAveraging!="element_box_phumob" ||
                   usesFermiDirac_ || !recombination_.srhEnabled() || recombination_.augerEnabled() || recombination_.bandToBandEnabled() || impactIonizationCoupled_ || electronQuantumPotentialConfig_.enabled)
                    throw std::runtime_error("Unqualified split operator physics");
                nlohmann::json j={{"schema","vela.split-dd-operator.v1"},{"temperature_K",300.},{"potential_scale",potentialScale}};
                for(const auto& a:pn)j["nodes"].push_back({a.psi,a.n,a.p,a.xpsi,a.xn,a.xp,a.potentialScale,a.eref,a.href,a.ni,a.vt,a.doping,a.voln,a.volp,a.vold,a.q,a.area,a.scale,a.interfaceRhs});
                for(const auto& e:pe)j["poisson_edges"].push_back({e.i,e.j,e.g});
                const auto& p=mobilityConfig_.phuMob;
                const auto& en=p.donorSpecies==PhuMobDonorSpecies::Arsenic?p.electronArsenic:p.electronPhosphorus;
                const auto& hp=p.holeBoron;
                j["phumob"]={{"electron",{en.muMax,en.muMin,en.theta,en.nRef,en.alpha}},{"hole",{hp.muMax,hp.muMin,hp.theta,hp.nRef,hp.alpha}},
                  {"common",{p.donorClusterReference,p.acceptorClusterReference,p.donorClusterCoefficient,p.acceptorClusterCoefficient,p.electronMassRatio,p.holeMassRatio,p.conwellWeisskopfFactor,p.brooksHerringFactor,p.electronHoleScatteringFactor,p.holeElectronScatteringFactor,p.gA,p.gB,p.gC,p.gAlpha,p.gAlphaPrime,p.gBeta,p.gGamma,p.internalConcentrationToCm3,p.internalMobilityToCm2PerVS}}};
                const Real cs=scaling_.enabled?scaling_.C0*scaling_.D0:1.;
                for(int i=0;i<N;++i) {
                    Real concentration=recombination_.srhDopingConcentration(doping_.donors(i),doping_.acceptors(i));
                    j["node_physics"].push_back({doping_.donors(i),doping_.acceptors(i),recombination_.electronLifetime(concentration),recombination_.holeLifetime(concentration),vol[i]*sourceIntegralFactor/cs});
                    j["baseline_density"].push_back({n(i),p(i)});
                }
                j["boundary"]=nlohmann::json::array({nlohmann::json::array(),nlohmann::json::array(),nlohmann::json::array()});
                for(const auto& [node,z]:bcs.psi)j["boundary"][0].push_back({node,z});
                for(const auto& [node,z]:bcs.phin)j["boundary"][1].push_back({node,z-electronQuasiFermiReferenceAt(node)/potentialScale});
                for(const auto& [node,z]:bcs.phip)j["boundary"][2].push_back({node,z-holeQuasiFermiReferenceAt(node)/potentialScale});
                const detail::EdgeMobilityCarrierState populations{0,0,0,0,&n,&p};
                for(Index e=0;e<mesh_.numEdges();++e) {
                    const auto& edge=mesh_.getEdge(e);if(edge.length<1e-30)continue;
                    nlohmann::json weights=nlohmann::json::array(),ports=nlohmann::json::object();
                    long double geometry=0;
                    for(Index cid:edgeCells_[e]) {const auto& m=cellMaterials_[cid];if(m.ni<=0 && m.mun<=0 && m.mup<=0)continue;geometry+=detail::cellBoxEdgeCoefficient(mesh_,cid,e);}
                    if(geometry>0)for(Index cid:edgeCells_[e]) {
                        const auto& m=cellMaterials_[cid];if(m.ni<=0 && m.mun<=0 && m.mup<=0)continue;
                        if(m.temperature_K.value_or(300.)!=300. || m.mun<=0 || m.mup<=0)throw std::runtime_error("Unqualified split mobility material");
                        Real g=detail::cellBoxEdgeCoefficient(mesh_,cid,e);if(g==0)continue;
                        const auto measures=detail::cellBoxNodeMeasures(mesh_,cid);long double volume=(long double)measures[0]+measures[1]+measures[2];
                        for(int k=0;k<3;++k)weights.push_back({mesh_.getCell(cid).node_ids[k],static_cast<Real>((long double)g*measures[k]/volume/geometry)});
                    }
                    if(couple[e]>0)for(const auto& c:mesh_.contacts()) {
                        bool i=std::find(c.node_ids.begin(),c.node_ids.end(),edge.n0)!=c.node_ids.end();
                        bool jn=std::find(c.node_ids.begin(),c.node_ids.end(),edge.n1)!=c.node_ids.end();
                        if(i!=jn)ports[c.name]=i?1.:-1.;
                    }
                    const Real mn=detail::elementBoxEdgeMobility(edgeCells_,mesh_,doping_,*mobility_,cellMaterials_,e,CarrierType::Electron,mobilityConfig_,&populations);
                    const Real mp=detail::elementBoxEdgeMobility(edgeCells_,mesh_,doping_,*mobility_,cellMaterials_,e,CarrierType::Hole,mobilityConfig_,&populations);
                    j["transport_edges"].push_back({{"id",e},{"i",edge.n0},{"j",edge.n1},{"coefficient",Vt_*fieldFactor*couple[e]/edge.length/cs},{"weights",weights},{"ports",ports},{"baseline_mu",{mn,mp}}});
                }
                j["current_factor"]=constants::q*cs*(scaling_.enabled?scaling_.currentDensityLineIntegralFactor:1.)*1e-6;
                std::ofstream(path)<<std::setprecision(17)<<j.dump(2);
            }
        }
'''
# avoid shadowing the hole population vector with the PhuMob parameter alias
HOOK=HOOK.replace('const auto& p=mobilityConfig_.phuMob;','const auto& pm=mobilityConfig_.phuMob;').replace('p.donor','pm.donor').replace('p.electron','pm.electron').replace('p.hole','pm.hole').replace('p.acceptor','pm.acceptor').replace('p.conwell','pm.conwell').replace('p.brooks','pm.brooks').replace('p.g','pm.g').replace('p.internal','pm.internal')

DRIVER=r'''
#include "vela/equation/SplitDDOperator.h"
#include <fstream>
#include <iomanip>
using namespace vela::split_dd;
int main(int argc,char**argv) {
    if(argc!=4)return 2;
    nlohmann::json data,jobs;std::ifstream(argv[1])>>data;std::ifstream(argv[2])>>jobs;
    Operator op(data);std::ofstream out(argv[3]);out<<std::setprecision(110)<<"case,kind,index,a,b,c,d,e,f,g,h,i\n";
    for(const auto& job:jobs) {
        State state=State::restore(job["state"],data["mesh_fingerprint"].get<std::string>());
        if(job.contains("step"))state=state.shifted(job["step"].get<std::vector<double>>(),job["alpha"]);
        if(job.value("repartition",false)) {
            auto j=state.checkpoint();for(auto& p:j["coordinates"]) {
                double hi=p[0],lo=p[1];double h=std::nextafter(hi,INFINITY);
                Wide low=Wide(hi)+Wide(lo)-Wide(h);double l=static_cast<double>(low);
                if(Wide(l)==low)p={h,l};
            }
            state=State::restore(j,data["mesh_fingerprint"].get<std::string>());
        }
        auto restored=State::restore(nlohmann::json::parse(state.checkpoint().dump()),data["mesh_fingerprint"].get<std::string>());
        for(std::size_t k=0;k<state.size()*3;++k)if(restored.coordinate(k)!=state.coordinate(k))return 3;
        auto r=op.evaluate(restored);std::string label=job["label"];
        for(std::size_t i=0;i<state.size();++i)out<<label<<",node,"<<i<<','<<r.n[i]<<','<<r.p[i]<<','<<r.mun[i]<<','<<r.mup[i]<<','<<r.srh[i]<<','<<r.residual[i]<<','<<r.residual[state.size()+i]<<','<<r.residual[2*state.size()+i]<<','<<r.psi[i]<<'\n';
        for(std::size_t i=0;i<r.electronFlux.size();++i)out<<label<<",edge,"<<i<<','<<r.electronFlux[i]<<','<<r.holeFlux[i]<<",,,,,,,\n";
        for(const auto& [name,z]:r.currents)out<<label<<",port,"<<name<<','<<z<<",,,,,,,,\n";
        std::ofstream(std::string(argv[3])+"."+label+".json")<<state.checkpoint().dump(2);
    }
}
'''

def build():
    import shutil
    original=LOCAL.parent/'probe'; old=read(OUT.parent/'build_evidence.json')['input_hashes']
    dest=LOCAL/'probe';dest.mkdir(parents=True,exist_ok=False)
    reused=[]
    for name in ('runner.exe','assembler.cpp','assembler.o','compile.json','compile.log','link.json','link.log'):
        p=original/name;assert prior.sha(p)==old[str(p.relative_to(ROOT)).replace('\\','/')]
        shutil.copy2(p,dest/name);reused.append(p)
    (dest/'driver.cpp').write_text(DRIVER);cmd=['D:/msys64/ucrt64/bin/c++.exe','-std=c++20','-O2','-I'+str(ROOT/'include'),str(dest/'driver.cpp'),'-o',str(dest/'driver.exe')]
    write(dest/'driver_compile.json',cmd);p=subprocess.run(cmd,env=env(),capture_output=True,text=True);(dest/'driver.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),OUT.parent/'build_evidence.json',ROOT/'include/vela/numerics/SplitDDState.h',ROOT/'include/vela/equation/SplitDDOperator.h',ROOT/'include/vela/numerics/SplitCoordinate.h']+reused+[p for p in dest.iterdir() if p.is_file()])

def run():
    verify(OUT/'build_evidence.json');verify(prior.OUT/'replay_evidence.json');records=[];files=[OUT/'build_evidence.json',prior.OUT/'replay_evidence.json']
    for c in rows(prior.OUT/'identity.csv'):
        base=Path(c['new_dest']);dest=LOCAL/'states'/c['device']/c['vd']/c['index'];dest.mkdir(parents=True,exist_ok=False)
        cfg=read(base/'config.json');cfg.update(residual_output_csv=str(dest/'baseline.csv'),contact_edge_output_csv=str(dest/'contact.csv'));write(dest/'config.json',cfg)
        e=env();e['VELA_CONTRACT_EXPORT']=str(dest/'operator_raw.json')
        p=subprocess.run([str(LOCAL/'probe/runner.exe'),'--config',str(dest/'config.json'),'--log','off'],cwd=ROOT,env=e,capture_output=True,text=True);(dest/'probe.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-1500:]
        assert prior.sha(dest/'baseline.csv')==prior.sha(base/'residual.csv')
        data=read(dest/'operator_raw.json');data['mesh_fingerprint']=prior.sha(Path(cfg['mesh_file']));write(dest/'operator.json',data)
        seed=rows(cfg['state_file']);n=len(seed);keys=('packed_psi','packed_electron_qf_increment','packed_hole_qf_increment')
        state=dict(schema='vela.split-dd-state.v1',mesh_fingerprint=data['mesh_fingerprint'],potential_scale_V=float(seed[0]['packed_potential_scale_V']),electron_reference_V=[float(r['electron_qf_reference_V']) for r in seed],hole_reference_V=[float(r['hole_qf_reference_V']) for r in seed],coordinates=[[float(r[k]),0.] for k in keys for r in seed])
        original=rows(Path(c['dest'])/'trial_0.csv');direction=[float(r['capped_step']) for r in original]
        jobs=[dict(label='base',state=state),dict(label='partition',state=state,repartition=True)]
        for alpha in (1.,2.**-12):
            for sign in (-1,1):jobs.append(dict(label=('plus' if sign>0 else 'minus')+('_full' if alpha==1 else '_small'),state=state,step=direction,alpha=sign*alpha))
        # One-ULP/64 probe in each block at the final worst-hole node.
        for b in range(3):
            step=[0.]*(3*n);node=int(read(Path(c['dest'])/'config.status.json')['carrier_row_convergence']['max_ratio_node']);step[b*n+node]=math.ulp(state['coordinates'][b*n+node][0])/64
            jobs.append(dict(label=f'micro_{b}',state=state,step=step,alpha=1.))
        write(dest/'jobs.json',jobs)
        p=subprocess.run([str(LOCAL/'probe/driver.exe'),str(dest/'operator.json'),str(dest/'jobs.json'),str(dest/'results.csv')],env=env(),capture_output=True,text=True);(dest/'driver.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-1500:]
        rec=dict(device=c['device'],vd=c['vd'],index=c['index'],dest=str(dest),old_dest=c['dest'],baseline_identity=True,jobs=len(jobs));records.append(rec);print(rec,flush=True)
        files += [Path(cfg['state_file'])]+[p for p in dest.iterdir() if p.is_file()]
    csvout(OUT/'runs.csv',records);freeze(OUT/'run_evidence.json',files+[OUT/'runs.csv'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
