#!/usr/bin/env python3
"""Audit ligand coverage and adjacency between predicted pocket voxel sets.

The ligand is used only for post-prediction evaluation. No merges are performed.
"""
import argparse
import json
import math
from collections import defaultdict


def read_atoms(path, ligand, chain, residue):
    found = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith("HETATM") or line[17:20].strip() != ligand:
            continue
        if chain is not None and line[21:22].strip() != chain:
            continue
        if residue is not None and line[22:26].strip() != str(residue):
            continue
        element = line[76:78].strip().upper()
        if element in ("H", "D"):
            continue
        found.append({"name": line[12:16].strip(), "chain": line[21:22].strip(),
                      "residue": line[22:26].strip(),
                      "xyz": tuple(float(line[i:i + 8]) for i in (30, 38, 46))})
    instances = {(a["chain"], a["residue"]) for a in found}
    if len(instances) != 1 or not found:
        raise ValueError("Ligand instance missing or ambiguous; specify chain and residue number")
    return found


def read_pockets(path):
    groups = defaultdict(list)
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith("HETATM") or line[17:20].strip() != "PKT":
            continue
        pid = int(line[22:26])
        groups[pid].append(tuple(float(line[i:i + 8]) for i in (30, 38, 46)))
    if not groups:
        raise ValueError("No pocket pseudoatoms in predicted PDB")
    return groups


def grid_adjacency(groups, pitch):
    # Exact integer lattice site for input pocket coordinates, robust to PDB rounding.
    origin = tuple(min(pt[d] for group in groups.values() for pt in group) for d in range(3))
    lattice = {}
    for pid, points in groups.items():
        for pt in points:
            key = tuple(round((pt[d] - origin[d]) / pitch) for d in range(3))
            if key in lattice and lattice[key] != pid:
                raise ValueError("Overlapping pocket voxel centers")
            lattice[key] = pid
    faces = defaultdict(int)
    for key, pid in lattice.items():
        for d in range(3):
            neigh = list(key)
            neigh[d] += 1
            other = lattice.get(tuple(neigh))
            if other is not None and other != pid:
                faces[tuple(sorted((pid, other)))] += 1
    return faces


def audit(ligand_atoms, pockets, distance, pitch):
    cutoff_sq = distance * distance
    atom_rows = []
    hit_by_pocket = defaultdict(set)
    for atom in ligand_atoms:
        hits = []
        nearest = None
        for pid, points in pockets.items():
            d2 = min(sum((a-b)**2 for a,b in zip(atom["xyz"], pt)) for pt in points)
            if nearest is None or d2 < nearest[1]:
                nearest = (pid, d2)
            if d2 <= cutoff_sq:
                hits.append({"pocket": pid, "distance_a": round(math.sqrt(d2), 4)})
                hit_by_pocket[pid].add(atom["name"])
        atom_rows.append({"atom": atom["name"], "nearest_pocket": nearest[0],
                          "nearest_distance_a": round(math.sqrt(nearest[1]), 4),
                          "covered_by": sorted(hits, key=lambda x: x["distance_a"])})
    faces = grid_adjacency(pockets, pitch)
    adjacency = [{"a": a, "b": b, "contact_faces": n,
                  "contact_area_a2": round(n*pitch*pitch, 3)}
                 for (a,b),n in faces.items()]
    adjacency.sort(key=lambda e: (-e["contact_faces"], e["a"], e["b"]))
    ranking = sorted(hit_by_pocket, key=lambda p: (-len(hit_by_pocket[p]), p))
    top = ranking[0] if ranking else None
    absent = [row for row in atom_rows if not any(v["pocket"] == top for v in row["covered_by"])] if top else atom_rows
    return {"ligand_heavy_atoms": len(ligand_atoms), "coverage_distance_a": distance,
            "best_coverage_pocket_id": top, "best_covered_atoms": len(hit_by_pocket[top]) if top else 0,
            "missing_from_best": absent, "atoms": atom_rows, "pocket_contact_graph": adjacency,
            "note": "Ligand used only for evaluation. Shared grid faces do not establish valid physical basin merging."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--ligand-chain")
    parser.add_argument("--ligand-resseq")
    parser.add_argument("--pockets-pdb", required=True)
    parser.add_argument("--distance", type=float, default=3.0)
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--output")
    args = parser.parse_args()
    if not math.isfinite(args.distance) or args.distance <= 0 or not math.isfinite(args.spacing) or args.spacing <= 0:
        parser.error("distance and spacing must be finite and positive")
    result = audit(read_atoms(args.crystal, args.ligand, args.ligand_chain, args.ligand_resseq),
                   read_pockets(args.pockets_pdb), args.distance, args.spacing)
    encoded = json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            stream.write(encoded + "\n")
    else:
        print(encoded)


if __name__ == "__main__":
    main()
