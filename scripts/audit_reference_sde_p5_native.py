"""Audit controlled P5 Sentaurus outputs without storing vendor files in Git."""
from pathlib import Path
import hashlib, json, math, re, sys
import argparse
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'))
from compare_sentaurus_tdr_tdx import parse_grd
from vela_sdevice.mesh_import import _triangle_from_edges



def grid(name):
    g=parse_grd(P/(name+'.grd'));local={};areas={}
    for region,ids in g['region_elements'].items():
        vertices=set();area=0.
        for i in ids:
            el=g['elements'][i]
            if el[0]!=2:continue
            ns=_triangle_from_edges(el,g['edges']);vertices.update(ns)
            a,b,c=[g['vertices'][j] for j in ns]
            area+=abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
        if area:local[region]=sorted(vertices);areas[region]=area
    return g,local,areas


def datasets(name):
    text=(P/(name+'.dat')).read_text()
    pat=r'Dataset\s*\("([^"]+)"\)\s*\{(.*?)Values\s*\((\d+)\)\s*\{(.*?)\n\s*\}'
    for m in re.finditer(pat,text,re.S):
        validity=re.search(r'validity\s*=\s*\[([^]]*)',m[2])[1]
        regions=re.findall(r'"([^"]+)"',validity)
        assert len(regions)==1,regions
        values=[float(v) for v in m[4].split()];assert len(values)==int(m[3])
        yield m[1],regions[0],values


def expected_bv(x,y,region):
    if region!='R.Substrate':return 0.,0.
    donor=0.
    for peak,depth,factor,start in [(5e20,.1,.4,.165),(1e19,.05,.8,.065)]:
        sigma=depth/math.sqrt(2*math.log(peak/1e18))
        donor+=peak*math.exp(-.5*(y/sigma)**2)*.5*(math.erf((abs(x)-start)/(factor*sigma))-math.erf((abs(x)-1.13)/(factor*sigma)))
    return donor,1e18


def main():
    global P
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    P=args.input_dir
    report={'schema':'vela.sde_p5_native_oracle.v1','vendor_version':'U-2022.12','status':'passed',
            'context':{'bv_source_modified':False,'wrapper_defines':{'+z':'+z'},'raw_source_execution':'rejected: unbound +z'},
            'native_topology_equivalence':False,'tests':{},'input_sha256':{}}
    for name in ['native_erf','native_removed','native_bv','native_ratio','native_min','native_max']:
        for ext in ['grd','dat']:
            p=P/(name+'.'+ext)
            report['input_sha256'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        g,local,areas=grid(name)
        record={'vertices':g['vertex_count'],'triangles':sum(e[0]==2 for e in g['elements']),'volume_area_um2':areas}
        if name in ['native_erf','native_bv']:
            maximum=0.;absolute=0.;count=0;seen=set()
            for field,region,values in datasets(name):
                seen.add((field,region))
                ids=local[region];assert len(ids)==len(values),(name,region,len(ids),len(values))
                for i,value in zip(ids,values):
                    x,y=g['vertices'][i]
                    if name=='native_erf':
                        sigma=.1/math.sqrt(2*math.log(100))
                        donor=1e18*math.exp(-.5*(y/sigma)**2)*.5*(math.erf(x/(.4*sigma))-math.erf((x-1)/(.4*sigma)));acceptor=0.
                        scale=1e18
                    else:donor,acceptor=expected_bv(x,y,region);scale=5e20
                    expected={'DopingConcentration':donor-acceptor,'ArsenicActiveConcentration':donor,'BoronActiveConcentration':acceptor}[field]
                    error=abs(value-expected)
                    # Report peak-normalized error; native erf cancellation in tiny
                    # tails is also recorded as an absolute error, never hidden.
                    maximum=max(maximum,error/scale);absolute=max(absolute,error);count+=1
            required = ({('DopingConcentration','R.Si'),('ArsenicActiveConcentration','R.Si')} if name=='native_erf'
                        else {(field,'R.Substrate') for field in ('DopingConcentration','ArsenicActiveConcentration','BoronActiveConcentration')})
            assert required <= seen and count>0, 'missing native oracle datasets'
            record.update(checked_region_local_values=count,maximum_peak_normalized_error=maximum,maximum_absolute_error_cm3=absolute)
            assert maximum<=1e-6,record
        if name=='native_removed':assert math.isclose(sum(areas.values()),.82,abs_tol=1e-12)
        report['tests'][name]=record
    one,low,high=(report['tests'][n]['vertices'] for n in ['native_ratio','native_min','native_max'])
    assert one==high and low>one
    report['ratio_one_rule']='maximum size (live controls)'
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'ratio_one_rule':report['ratio_one_rule']}))


if __name__=='__main__':
    main()
