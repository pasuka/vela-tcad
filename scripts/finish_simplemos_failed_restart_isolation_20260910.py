"""Finish the probe interrupted by Windows path length; preserve its first output."""
from pathlib import Path
import isolate_simplemos_production_failed_restart_20260910 as old
v,a,d,run=old.v,old.a,old.d,old.run

def main():
    v.configure();a.verify(old.OUT/'freeze.json')
    jobs=a.read(old.OUT/'contract.json')['jobs'];job=jobs[-1]
    assert job['model']=='unprojected_contacts'
    run.LOCAL=v.REPO/'build-release/pp_iso';run.OUT=v.OUT
    row=run.attempt(job,Path(job['seed']),0)
    rows=[a.read(old.LOCAL/'dc'/j['model']/j['case']/'diagnostic/vg_000/attempt_0/result.json') for j in jobs[:-1]]+[row]
    run.v.csv_union(old.OUT/'results.csv',rows)
    a.write(old.OUT/'completion.json',dict(reason='Last acceptance probe produced CSV but its long stdout filename exceeded Windows path support. Rerun only this identical diagnostic control under build-release/pp_iso; original output retained.',cohort_credit=False))
    d.matrix.freeze(old.OUT/'evidence.json',[Path(__file__).resolve(),old.OUT/'freeze.json',old.OUT/'completion.json',old.OUT/'results.csv']+[p for root in (old.LOCAL,run.LOCAL) for p in root.rglob('*') if p.is_file()])

if __name__=='__main__':main()
