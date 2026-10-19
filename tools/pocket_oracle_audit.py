#!/usr/bin/env python3
"""Crystal-only diagnostic: distinguish discovery, segmentation and ranking failures."""
import argparse
import json
from pathlib import Path
from pocket_contact_audit import read_atoms, read_pockets

def audit(pdb, crystal, chain, seq, cutoff=3.0):
    atoms=read_atoms(str(crystal),"STI",chain,seq)
    pockets=read_pockets(str(pdb))
    hits={}
    for idx,coords in pockets.items():
        hits[idx]={i for i,atom in enumerate(atoms)
            if any(sum((atom["xyz"][k]-p[k])**2 for k in range(3))<=cutoff**2
                   for p in coords)}
    covered=set().union(*hits.values()) if hits else set()
    remaining=set(covered)
    selection=[]
    while remaining:
        best=max(hits,key=lambda k:(len(hits[k]&remaining),-len(pockets[k]),-k))
        newly=hits[best]&remaining
        if not newly:break
        selection.append({"pocket":best,"new_atoms":len(newly),
                          "covered_atoms_total":len(covered)-len(remaining)+len(newly)})
        remaining-=newly
    return {"chain":chain,"residue":seq,"pocket_count":len(pockets),
            "ligand_heavy_atoms":len(atoms),"oracle_union_covered":len(covered),
            "oracle_union_coverage":len(covered)/len(atoms),
            "uncovered_atom_names":[atom["name"] for i,atom in enumerate(atoms)
                                     if i not in covered],
            "greedy_set_cover":selection,
            "top_pocket_contributions":[{"pocket":k,"atoms_covered":len(v),
                                        "voxel_count":len(pockets[k])}
                for k,v in sorted(hits.items(),key=lambda kv:-len(kv[1]))[:15]]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crystal",type=Path,required=True)
    p.add_argument("--pockets-a",type=Path,required=True)
    p.add_argument("--pockets-b",type=Path,required=True)
    p.add_argument("--output",type=Path,default=Path("1iep_oracle_audit.json"))
    p.add_argument("--distance",type=float,default=3.0)
    args=p.parse_args()
    data={"note":"Retrospective oracle diagnostics only; not valid as prediction input",
        "A":audit(args.pockets_a,args.crystal,"A","201",args.distance),
        "B":audit(args.pockets_b,args.crystal,"B","202",args.distance)}
    args.output.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
if __name__=="__main__":main()
