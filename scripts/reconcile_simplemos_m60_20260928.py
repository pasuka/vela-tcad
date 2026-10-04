"""Read-only join of the frozen September 27 Vela run and September 2 M60.

No simulator is invoked and neither reference nor acceptance policy is changed.
"""
import csv, hashlib, json, math, re, shutil, statistics
from pathlib import Path
from sentaurus_import import parse_quoted_list, parse_values_block
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'reference_tcad/simplemos_sentaurus2022'
B=ROOT/'build/outlier_analysis_20260928'
M=R/'tight_convergence_port_burst'
RAW=ROOT/'build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst'
O=ROOT/'build/outlier_m60_review_20260928'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(name,data):
    with (O/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in data for k in r)));w.writeheader();w.writerows(data)
def key(r):
    vg=float(r.get('vg',r.get('gate_voltage_V')));idx=round(vg/.05)
    assert abs(vg-idx*.05)<1e-10
    return r['case'],idx
def keyed(data):
    result={key(r):r for r in data};assert len(result)==len(data);return result
def remove_blocks(text,name):
    while m:=re.search(r'\b'+name+r'\s*\(',text):
        level=1;i=m.end()
        while level:
            level+=(text[i]=='(')-(text[i]==')');i+=1
        text=text[:m.start()]+text[i:]
    return text
def blocks(text,name):
    out=[]
    for m in re.finditer(r'\b'+name+r'(?:\([^)]*\))?\s*\{',text):
        level=1;i=m.end()
        while level:
            level+=(text[i]=='{')-(text[i]=='}');i+=1
        out.append(text[m.start():i])
    return out
def canonical(text):
    text=re.sub(r'(?<![\w])[-+]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[eE][-+]?\d+)?',lambda m:format(float(m[0]),'.12g'),text)
    return re.sub(r'\s+','',text)
def main():
    O.mkdir(exist_ok=True);verified=[]
    ev=read(R/'simplemos_m60_tight_convergence_port_burst_evidence.json')
    for section in ['source_hashes','implementation_hashes','artifacts']:
        for name,expected in ev[section].items():
            p=ROOT/name;assert p.exists(),name;assert sha(p)==expected,name
            verified.append(dict(path=name,sha256=expected,group=section))
    for name,expected in read(B/'cloud_hashes.json').items():assert sha(B/name)==expected,name
    v=keyed(rows(B/'all_points.csv'));p=keyed(rows(M/'m60_tight_default_direct_point_ledger.csv'))
    assert len(v)==816 and v.keys()==p.keys()
    terminal=rows(M/'m60_tight_terminal_component_ledger.csv')
    tt={(key(r),r['algorithm'],r['contact'],r['component']):float(r['current_A_per_um']) for r in terminal}
    assert len(tt)==len(terminal)==816*2*4*3
    def current(k,alg,c,comp):return tt[k,alg,c,comp]
    manifests=read(RAW/'sentaurus_manifest.json')['cases'];logs={r['run_case']:r for r in rows(M/'m60_cnormprint_log_ledger.csv')}
    raw_checks=[];deck_checks=[]
    for c in manifests:
        deck=ROOT/c['deck'];tdr=ROOT/c['input_tdr']
        assert sha(deck)==c['deck_sha256'];assert sha(tdr)==c['input_tdr_sha256']
        old=ROOT/f"build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle/{c['device']}/input_fps.tdr"
        assert sha(tdr)==sha(old)
        raw=RAW/'sentaurus_raw/sentaurus_bundle'/c['device']
        log=raw/c['expected_console'];assert sha(log)==logs[c['run_case']]['console_sha256']
        assert logs[c['run_case']]['completed_without_exit_failure']=='1'
        text=(raw/c['expected_current']).read_text();names=parse_quoted_list(text,'datasets')
        rawrows=[dict(zip(names,r)) for r in parse_values_block(text,len(names))]
        for k in [k for k in p if k[0]==c['base_case']]:
            match=[r for r in rawrows if abs(r['gate OuterVoltage']-k[1]*.05)<1e-10];assert len(match)==1
            for contact in ['source','drain','substrate','gate']:
                for comp,label in [('electron','eCurrent'),('hole','hCurrent'),('total','TotalCurrent')]:
                    assert match[0][contact+' '+label]==current(k,c['algorithm'],contact,comp)
        raw_checks.append(dict(run_case=c['run_case'],plt_sha256=sha(raw/c['expected_current']),console_sha256=sha(log),points=51))
        baseline=ROOT/f"build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle/{c['device']}/{c['base_case']}_des.cmd"
        a=baseline.read_text();b=deck.read_text()
        physics=canonical(''.join(blocks(a,'Physics')))==canonical(''.join(blocks(b,'Physics')))
        contacts=canonical(''.join(blocks(a,'Electrode')))==canonical(''.join(blocks(b,'Electrode')))
        solve=canonical(remove_blocks(''.join(blocks(a,'Solve')),'Plot'))==canonical(remove_blocks(''.join(blocks(b,'Solve')),'Plot'))
        assert physics and contacts and solve,c['run_case']
        deck_checks.append(dict(run_case=c['run_case'],physics_equal=physics,contacts_equal=contacts,solve_equal_except_Plot=solve,tdr_equal=True,math=blocks(b,'Math')[0]))
    joined=[];selected=[];vela_ports={(key(r),r['contact']):r for r in rows(B/'vela_terminals.csv')}
    for k,r in v.items():
        s=p[k];vd=float(r['vd']);assert float(r['native_Id_A_per_um'])==float(s['m46_default_drain_current_A_per_um'])
        new=float(r['current_A_per_um']);loose=float(r['native_Id_A_per_um']);tight=float(s['tight_default_drain_current_A_per_um'])
        assert tight==current(k,'default','drain','total')
        assert float(s['tight_direct_drain_current_A_per_um'])==current(k,'direct','drain','total')
        kcl=math.fsum(current(k,'default',c,'total') for c in ['source','drain','substrate','gate'])
        gap=current(k,'default','substrate','electron')-current(k,'direct','substrate','electron')
        row=dict(case=k[0],index=k[1],device=r['device'],vd=vd,vg=float(r['vg']),
            vela_Id_A_per_um=new,original_Id_A_per_um=loose,tight_default_Id_A_per_um=tight,
            original_error_percent=100*(new/loose-1),tight_error_percent=100*(new/tight-1),
            tight_native_kcl_A_per_um=kcl,tight_native_kcl_percent=100*kcl/abs(tight),
            tight_substrate_observer_gap_A_per_um=gap,tight_substrate_observer_gap_percent=100*gap/abs(tight),
            proposed_1e_3_kcl_pass=abs(kcl/tight)<1e-3,proposed_1e_3_observer_pass=abs(gap/tight)<1e-3,
            sub_fA=abs(tight)<1e-15,vela_qualified=r['qualified'])
        joined.append(row)
        if abs(row['original_error_percent'])>2:
            vi=vela_ports[k,'substrate'];direct=current(k,'direct','substrate','electron');ratio=float(vi['electron_A_per_um'])/direct
            selected.append(dict(row,vela_sub_e_A_per_um=vi['electron_A_per_um'],tight_default_sub_e_A_per_um=current(k,'default','substrate','electron'),tight_direct_sub_e_A_per_um=direct,
                vela_over_tight_direct_sub_e=ratio,charge_normalized_ratio_if_manual_q=ratio/(1.602176634/1.602192),tight_direct_Id_A_per_um=current(k,'direct','drain','total')))
    table('joined_816.csv',joined);table('six_original_outliers.csv',selected)
    table('raw_plt_log_verification.csv',raw_checks);table('deck_equivalence.csv',deck_checks)
    table('verified_m60_evidence_hashes.csv',verified)
    # Extract actual accepted Newton histories at the six targets, not estimated iterations.
    stopping=[]
    for r in selected:
        log=RAW/f"sentaurus_raw/sentaurus_bundle/{r['device']}/m60_default_{r['case']}.console.log"
        text=log.read_text();text=text[text.index('Contact gate : 2.5V'):]
        starts=list(re.finditer(r'Computing step from t=([^ ]+) to t=([^ ]+) \(Stepsize: ([^)]+)\)',text))
        target=r['vg']/2.5
        match=[(i,m) for i,m in enumerate(starts) if abs(float(m[2])-target)<1e-10]
        assert len(match)==1
        i,m=match[0];block=text[m.start():starts[i+1].start() if i+1<len(starts) else len(text)]
        nr=[line.split() for line in block.splitlines() if re.match(r'^\s*\d+\s+\d\.\d+e[+-]\d+',line,re.I)]
        last=nr[-1];assert len(last)>=8
        reason=re.search(r'Finished, because\.\.\.\s*\n([^\n]+)',block)[1]
        name=r['case']+'_vg_'+format(r['vg'],'.2f')+'_newton.txt'
        (O/name).write_text(block,encoding='utf-8')
        stopping.append(dict(case=r['case'],vg=r['vg'],iterations=int(last[0]),rhs=float(last[1]),update_error=float(last[4]),stop_reason=reason.strip(),excerpt=name))
    table('six_strict_newton_stops.csv',stopping)
    fields=rows(M/'m60_default_direct_field_invariance_ledger.csv');states=sorted({r['state'] for r in fields})
    assert all(r['within_state_invariance']=='True' and float(r['maximum_absolute_difference'])==0 for r in fields)
    stats=dict(points=len(joined),exact_original_reference_matches=len(joined),old_outliers=len(selected),
        tight_outliers=sum(abs(r['tight_error_percent'])>2 for r in joined),
        max_abs_tight_error_percent=max(abs(r['tight_error_percent']) for r in joined),
        max_abs_tight_error_point=max(joined,key=lambda r:abs(r['tight_error_percent'])),
        high_vd_max_abs_tight_error_percent=max(abs(r['tight_error_percent']) for r in joined if r['vd']==1),
        low_vd_mean_signed_tight_error_percent=statistics.mean(r['tight_error_percent'] for r in joined if r['vd']==.05),
        proposed_rule_all_point_failure_count=sum(not(r['proposed_1e_3_kcl_pass'] and r['proposed_1e_3_observer_pass']) for r in joined),
        proposed_rule_sub_fA_failure_count=sum(r['sub_fA'] and not(r['proposed_1e_3_kcl_pass'] and r['proposed_1e_3_observer_pass']) for r in joined),
        proposed_rule_sub_fA_point_count=sum(r['sub_fA'] for r in joined),
        original_outliers_failing_proposed_rule=[r['device'] for r in selected if not(r['proposed_1e_3_kcl_pass'] and r['proposed_1e_3_observer_pass'])],
        observer_invariance_rows=len(fields),observer_invariance_states=states,
        m60_verified_hash_count=len(verified),raw_plts_verified=len(raw_checks),native_component_values_verified=len(terminal),
        raw_logs_verified=len(raw_checks),deck_pairs_verified=len(deck_checks),solver_invoked=False,reference_replaced=False,policy_changed=False)
    (O/'summary.json').write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(stats,indent=2,ensure_ascii=False))
    for r in selected:print(r['device'],r['tight_error_percent'],r['tight_native_kcl_percent'],r['tight_substrate_observer_gap_percent'],r['charge_normalized_ratio_if_manual_q'])

if __name__=='__main__':main()
