"""Replay low-field element mobility on an exported native temperature field.

Diagnostic only: native densities, element electric field and orientation are
inputs. No Newton solve, fit, or electrical/thermal qualification is performed.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from collections import defaultdict
import numpy as np
from audit_templates_ldmos_g3_idvg_shift_kcl import scalar_field, distribution
from audit_templates_ldmos_averagebox_node4492 import parse_debug_block


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("mesh", "export", "measure-debug", "parameter-audit", "probe", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--high-field-node-export",type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    mesh=read(args.mesh)
    xy={n["id"]:np.array([n["x"],n["y"]]) for n in mesh["nodes"]}
    regions={r["id"]:r["material"] for r in mesh["regions"]}
    cells=[c for c in mesh["triangles"] if regions[c["region_id"]]=="Si"]
    # Explicit topology and coordinate identity, rather than matching row counts.
    with (args.export/"nodes.csv").open() as f:
        native_nodes=list(csv.DictReader(f))
    for row in native_nodes:
        i=int(row["id"])
        if max(abs(xy[i][0]-float(row["x_um"])),abs(xy[i][1]-float(row["y_um"])))>1e-10:
            raise ValueError("Native and mesh coordinates differ")
    edges=defaultdict(list)
    for c in mesh["triangles"]:
        ns=c["node_ids"]
        for k in range(3): edges[tuple(sorted((ns[k],ns[(k+1)%3])))].append(c)
    interfaces=[(xy[a],xy[b]-xy[a]) for (a,b),adj in edges.items()
                if {regions[c["region_id"]] for c in adj}=={"Si","SiO2"}]
    starts=np.array([a for a,b in interfaces]); vec=np.array([b for a,b in interfaces])
    lengths=np.sum(vec*vec,axis=1)
    names=("DonorConcentration","AcceptorConcentration","eDensity","hDensity","LatticeTemperature","NearestInterfaceOrientation","eQuasiFermiPotential","hQuasiFermiPotential")
    fields={n:scalar_field(args.export,n) for n in names}
    nodes={n for c in cells for n in c["node_ids"]}
    for name in names:
        if not nodes.issubset(fields[name]): raise ValueError("Incomplete field "+name)
    contact_nodes={n for c in mesh["contacts"] if c["name"] in ("source","drain","substrate") for n in c["node_ids"]}
    contact_edges={tuple(sorted(e)) for c in mesh["contacts"] for e in c["edge_node_ids"]}
    tangents={}
    for edge,adj in edges.items():
        if (len(adj)==1 and edge not in contact_edges) or any(regions[c["region_id"]]!="Si" for c in adj):
            a,b=[xy[n] for n in edge]; tangent=(b-a)/np.linalg.norm(b-a)
            for c in adj:
                if regions[c["region_id"]]=="Si": tangents[c["id"]]=tangent
    if args.high_field_node_export:
        # Native control uses the same loaded state; a closed reference differs
        # by at most the native final update roundoff, checked explicitly here.
        for name in ("LatticeTemperature","eQuasiFermiPotential","hQuasiFermiPotential","eDensity","hDensity"):
            control=scalar_field(args.high_field_node_export,name)
            for n in nodes:
                tolerance=1e-9*max(1.,abs(fields[name][n]))
                if abs(control[n]-fields[name][n])>tolerance: raise ValueError("Native HFS state differs: "+name)
    distances={}
    for n in nodes:
        delta=xy[n]-starts; t=np.clip(np.sum(delta*vec,axis=1)/lengths,0,1)
        distances[n]=float(np.min(np.linalg.norm(delta-t[:,None]*vec,axis=1)))
    with (args.export/"fields/ElectricField_region0_cells.csv").open() as f:
        ef={int(r["cell_id"]):np.array([float(r["component0"]),float(r["component1"])]) for r in csv.DictReader(f)}
    measures=parse_debug_block(args.measure_debug,"Measure")
    groups=read(args.parameter_audit)["parameter_groups"]
    env=dict(os.environ)
    if os.name=="nt": env["PATH"]=r"D:\msys64\ucrt64\bin"+os.pathsep+env.get("PATH","")
    results={}
    for index,carrier in enumerate(("electron","hole")):
        prefix=("e","h")[index]
        with (args.export/("fields/"+prefix+"Mobility_region0_cells.csv")).open() as f:
            expected={int(r["cell_id"]):float(r["component0"]) for r in csv.DictReader(f)}
        for treatment in (("vertex_temperature",) if args.high_field_node_export else ("vertex_temperature","cell_mean_temperature","frozen_300K")):
            batches=defaultdict(list)
            for cell in cells:
                ids=cell["node_ids"]; pts=np.array([xy[n] for n in ids])
                d=np.array([distances[n] for n in ids])
                gradient=np.linalg.solve(pts[1:]-pts[0],d[1:]-d[0])
                en=abs(float(np.dot(ef[cell["id"]],gradient)))*100.
                mean_t=sum(fields["LatticeTemperature"][n] for n in ids)/3.
                q=np.array([fields[prefix+"QuasiFermiPotential"][n] for n in ids])
                qfg=np.linalg.solve(pts[1:]-pts[0],q[1:]-q[0])*1e6
                drive=abs(float(np.dot(qfg,tangents[cell["id"]]))) if cell["id"] in tangents else np.linalg.norm(qfg)
                if any(n in contact_nodes for n in ids): drive=np.linalg.norm(ef[cell["id"]])*100.
                parallel=np.linalg.norm(ef[cell["id"]]-float(np.dot(ef[cell["id"]],gradient))*gradient)*100.
                for k,n in enumerate(ids):
                    family=str(int(fields["NearestInterfaceOrientation"][n]))
                    temperature=fields["LatticeTemperature"][n] if treatment=="vertex_temperature" else mean_t if treatment=="cell_mean_temperature" else 300.
                    batches[family].append(dict(id=3*cell["id"]+k,donors_m3=fields["DonorConcentration"][n]*1e6,
                        acceptors_m3=fields["AcceptorConcentration"][n]*1e6,electrons_m3=fields["eDensity"][n]*1e6,
                        holes_m3=fields["hDensity"][n]*1e6,normalField_V_per_m=en,
                        interfaceDistance_m=distances[n]*1e-6,temperature_K=temperature))
                    if args.high_field_node_export:
                        fraction=fields[prefix+"Density"][n]/(fields[prefix+"Density"][n]+1e12)
                        batches[family][-1]["drivingField_V_per_m"]=float(fraction*drive+(1.-fraction)*parallel)
            values={}; high_values={}
            for family,states in batches.items():
                group=groups[family]
                params={k:v[index] if len(v)==2 else v[0] for k,v in group.items()
                        if k not in ("EnormMinimum","me_over_m0","mh_over_m0")}
                params["mass"]=group[("me_over_m0","mh_over_m0")[index]][0]
                params["other_mass"]=group[("mh_over_m0","me_over_m0")[index]][0]
                stem=prefix+"_"+treatment+"_"+family
                source=args.output/(stem+"_input.json"); destination=args.output/(stem+"_output.json")
                cfg=dict(temperature_K=300.,carrier=carrier,parameters_cm=params,states_SI=states)
                if args.high_field_node_export:
                    cfg["high_field"]=dict(vsat300_m_per_s=(1.07e5,8.37e4)[index],beta300=(1.109,1.213)[index],
                        vsat_exponent=(.87,.52)[index],beta_exponent=(.66,.17)[index])
                write(source,cfg)
                with destination.open("w",encoding="utf-8") as out:
                    subprocess.run([str(args.probe.resolve()),str(source.resolve())],env=env,stdout=out,check=True)
                values.update({r["id"]:r["mobility_m2_per_Vs"]*1e4 for r in read(destination)["results"]})
                if args.high_field_node_export:
                    high_values.update({r["id"]:r["high_field_mobility_m2_per_Vs"]*1e4 for r in read(destination)["results"]})
            records=[]; high_numerator=defaultdict(float); high_denominator=defaultdict(float)
            for cell in cells:
                i=cell["id"]; raw=measures[i]; weights=[raw[0],raw[2],raw[1]]
                predicted=sum(values[3*i+k]*w for k,w in enumerate(weights))/sum(weights)
                if args.high_field_node_export:
                    high=sum(high_values[3*i+k]*w for k,w in enumerate(weights))/sum(weights)
                    a,b,c=[xy[n] for n in cell["node_ids"]]; area=abs(float(np.linalg.det(np.array([b-a,c-a]))))/2.
                    for n in cell["node_ids"]: high_numerator[n]+=high/area; high_denominator[n]+=1./area
                records.append(dict(cell_id=i,native_cm2_per_Vs=expected[i],predicted_cm2_per_Vs=predicted,
                                    relative_error=predicted/expected[i]-1.))
            with (args.output/(prefix+"_"+treatment+".csv")).open("w",newline="",encoding="utf-8") as f:
                writer=csv.DictWriter(f,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
            results[prefix+"_"+treatment]=distribution(r["relative_error"] for r in records)
            if args.high_field_node_export:
                control=scalar_field(args.high_field_node_export,prefix+"Mobility")
                rows=[dict(node_id=n,native_cm2_per_Vs=control[n],predicted_cm2_per_Vs=high_numerator[n]/high_denominator[n],
                           relative_error=high_numerator[n]/high_denominator[n]/control[n]-1.) for n in sorted(nodes)]
                with (args.output/(prefix+"_high_field_nodes.csv")).open("w",newline="",encoding="utf-8") as f:
                    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
                results[prefix+"_high_field_nodes"]=distribution(r["relative_error"] for r in rows)
    provenance={str(p):sha(p) for p in [args.mesh,args.measure_debug,args.parameter_audit,args.probe,Path(__file__),args.export/"field_manifest.json",args.export/"nodes.csv"]}
    provenance.update({str(p):sha(p) for p in (args.export/"fields").glob("*.csv")})
    if args.high_field_node_export:
        provenance.update({str(p):sha(p) for p in (args.high_field_node_export/"fields").glob("*.csv")})
    write(args.output/"summary.json",dict(scope="Prescribed native fields; local IALMob, optional alpha=0 hot HFS and inverse-area nodal projection. No electrical or D0 qualification.",
        silicon_cells=len(cells),silicon_nodes=len(nodes),temperature_range_K=[min(fields["LatticeTemperature"].values()),max(fields["LatticeTemperature"].values())],
        results=results,sha256=provenance))
    print(json.dumps(results,indent=2))


if __name__=="__main__":
    main()
