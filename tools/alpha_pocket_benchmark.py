#!/usr/bin/env python3
"""Reference-only evaluation of ligint-pocket-alpha results (not used in detection)."""
import argparse
import csv
import math
from pathlib import Path

def crystal_ligand(path, resname, chain, resseq):
    points = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("HETATM", "ATOM  ")) or len(line) < 54:
            continue
        if line[17:20].strip().upper() != resname.upper():
            continue
        if chain and line[21:22].strip() != chain:
            continue
        if resseq is not None and line[22:26].strip() != str(resseq):
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else line[12:16].strip()[0].upper()
        if element in ("H", "D"):
            continue
        try:
            points.append(tuple(float(line[k:k+8]) for k in (30,38,46)))
        except ValueError:
            continue
    if not points:
        raise ValueError("No heavy atoms found for the specified ligand residue")
    return points

def distance(a, b):
    return math.dist(a, b)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--ligand", required=True, help="PDB ligand residue name, e.g. MK1")
    parser.add_argument("--ligand-chain", default="")
    parser.add_argument("--ligand-resseq", type=int)
    parser.add_argument("--pockets", required=True, help="Alpha pocket summary CSV")
    parser.add_argument("--spheres", required=True, help="Alpha sphere pseudoatom PDB")
    parser.add_argument("--output", required=True)
    parser.add_argument("--hit-distance", type=float, default=4.0)
    args=parser.parse_args()
    ligand=crystal_ligand(args.crystal,args.ligand,args.ligand_chain,args.ligand_resseq)
    pocket_spheres={}
    for line in Path(args.spheres).read_text().splitlines():
        if not line.startswith("HETATM") or len(line)<54:
            continue
        try:
            rank=int(line[22:26])
            center=tuple(float(line[k:k+8]) for k in (30,38,46))
            pocket_spheres.setdefault(rank,[]).append(center)
        except ValueError:
            continue
    with open(args.pockets,newline="") as f:
        pockets=list(csv.DictReader(f))
    records=[]
    for item in pockets:
        rank=int(item["rank"])
        center=tuple(float(item[k]) for k in ("center_x","center_y","center_z"))
        center_distance=min(distance(center,p) for p in ligand)
        sphere_distance=min((distance(s,p) for s in pocket_spheres.get(rank,[]) for p in ligand),default=float("nan"))
        records.append({"rank":rank,"sphere_count":item["sphere_count"],
            "center_min_ligand_A":round(center_distance,4),
            "sphere_min_ligand_A":round(sphere_distance,4),
            "sphere_hit":int(math.isfinite(sphere_distance) and sphere_distance<=args.hit_distance)})
    with open(args.output,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=records[0].keys() if records else
            ["rank","sphere_count","center_min_ligand_A","sphere_min_ligand_A","sphere_hit"])
        writer.writeheader();writer.writerows(records)
    print(f"Ligand heavy atoms: {len(ligand)}, macropockets: {len(records)}")
    for k in (1,3,5,10):
        subset=records[:k]
        hit=any(r["sphere_hit"] for r in subset)
        print(f"Top-{k} ligand-site hit (sphere <= {args.hit_distance:g} Å): {int(hit)}")
    if records:
        best=min(records,key=lambda r: r["sphere_min_ligand_A"] if math.isfinite(r["sphere_min_ligand_A"]) else float("inf"))
        print(f"Closest sphere: pocket {best['rank']}, {best['sphere_min_ligand_A']} Å")

if __name__=="__main__":
    main()
