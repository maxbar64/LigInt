#!/usr/bin/env python3
"""Evaluate predicted pocket unions against a crystal ligand; never used for ranking."""
import argparse
import json
import math
from pocket_contact_audit import read_atoms, read_pockets

def benchmark(atoms, pockets, candidates, cutoff):
    coords=[a["xyz"] for a in atoms]
    rows=[]
    for rank,g in enumerate(candidates,1):
        members=g["members"]
        if not members or any(m not in pockets for m in members):
            raise ValueError("Group references absent pocket IDs")
        voxels=[point for member in members for point in pockets[member]]
        hits=sum(any((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2<=cutoff**2
                      for y in voxels) for x in coords)
        rows.append({"rank":rank,"members":members,"volume_a3":g.get("volume_a3"),
                     "score":g.get("score"),"covered_atoms":hits,
                     "ligand_heavy_atoms":len(coords),
                     "coverage":hits/len(coords)})
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crystal",required=True)
    p.add_argument("--ligand",required=True)
    p.add_argument("--ligand-chain")
    p.add_argument("--ligand-resseq")
    p.add_argument("--pockets-pdb",required=True)
    p.add_argument("--groups-json",required=True)
    p.add_argument("--distance",type=float,default=3.0)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    if not math.isfinite(a.distance) or a.distance<=0:p.error("Invalid distance")
    atoms=read_atoms(a.crystal,a.ligand,a.ligand_chain,a.ligand_resseq)
    pockets=read_pockets(a.pockets_pdb)
    with open(a.groups_json,encoding="utf-8") as f: report=json.load(f)
    if report.get("pocket_count")!=len(pockets):
        raise ValueError("Pocket count mismatch: regenerate all outputs together")
    rows=benchmark(atoms,pockets,report["groups"],a.distance)
    result={"distance_a":a.distance,"ligand_heavy_atoms":len(atoms),
            "candidate_count":len(rows),"best_coverage":max((r["coverage"] for r in rows),default=0),
            "first_full_coverage_rank":next((r["rank"] for r in rows if r["covered_atoms"]==len(atoms)),None),
            "results":rows}
    with open(a.output,"w",encoding="utf-8") as f:json.dump(result,f,indent=2);f.write("\n")
if __name__=="__main__":main()
