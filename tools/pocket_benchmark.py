#!/usr/bin/env python3
"""LigInt pocket benchmark: compare predicted PDB pseudoatoms with a crystal ligand.

Standard library only. No ligand positions are supplied to pocket discovery.
"""
import argparse
import json
import math
from pathlib import Path


def read_pdb(path, ligand_resname=None, ligand_chain=None, ligand_resseq=None, predictions=False):
    atoms = []
    with open(path, encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line[:6].strip() != "HETATM":
                continue
            if len(line) < 54:
                continue
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            name = line[12:16].strip()
            if element in ("H", "D") or (not element and name.startswith(("H", "D"))):
                continue
            if predictions:
                if resname != "PKT":
                    continue
            elif resname != ligand_resname or (ligand_chain is not None and chain != ligand_chain) or (
                ligand_resseq is not None and resseq != str(ligand_resseq)
            ):
                continue
            try:
                xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            except ValueError:
                continue
            if all(math.isfinite(x) for x in xyz):
                atoms.append((resseq if predictions else (chain, resseq), xyz))
    return atoms


def coverage(ligand, points, radius):
    if not ligand or not points:
        return 0.0
    r2 = radius * radius
    found = 0
    for _, xyz in ligand:
        if any(sum((a-b)**2 for a, b in zip(xyz, pos)) <= r2 for pos in points):
            found += 1
    return found / len(ligand)


def evaluate(crystal, ligand, pocket_pdb, pocket_json, ligand_chain, ligand_resseq, threshold):
    ligand_atoms = read_pdb(crystal, ligand, ligand_chain, ligand_resseq)
    if not ligand_atoms:
        raise ValueError("Ligand not found; specify its residue name, chain and sequence explicitly")
    # Fail closed if one residue name matches multiple ligand copies.
    ligand_ids = {key for key, _ in ligand_atoms}
    if len(ligand_ids) != 1:
        raise ValueError("Ambiguous ligand: specify --ligand-chain and --ligand-resseq")
    predicted = read_pdb(pocket_pdb, predictions=True)
    grouped = {}
    for pid, pos in predicted:
        grouped.setdefault(int(pid), []).append(pos)
    with open(pocket_json, encoding="utf-8") as stream:
        meta = json.load(stream)
    pockets = meta["pockets"]
    if sum(int(p["points"]) for p in pockets) != len(predicted):
        raise ValueError("Mismatch between predicted PDB and JSON pocket sizes")
    rows = []
    for p in pockets:
        pid = int(p["id"])
        points = grouped.get(pid, [])
        if len(points) != int(p["points"]):
            raise ValueError(f"PDB/JSON mismatch for pocket {pid}")
        hit2 = coverage(ligand_atoms, points, 2.0)
        hit3 = coverage(ligand_atoms, points, threshold)
        rows.append({
            "rank": pid, "coverage_2A": round(hit2, 6),
            "coverage_threshold": round(hit3, 6),
            "volume_a3": p["volume_a3"],
            "heuristic_score": p.get("heuristic_score", ""),
            "solvent_accessible_fraction": p.get("solvent_accessible_fraction", ""),
            "accessibility_probe_0_8": p.get("accessibility_probe_0_8", ""),
            "accessibility_probe_2_0": p.get("accessibility_probe_2_0", "")
        })
    best = max(rows, key=lambda r: (r["coverage_threshold"], -r["rank"]), default=None)
    return {
        "structure": str(crystal), "ligand": ligand,
        "ligand_id": next(iter(ligand_ids)),
        "ligand_heavy_atoms": len(ligand_atoms),
        "coverage_radius_a": threshold,
        "top1_coverage": rows[0]["coverage_threshold"] if rows else 0,
        "best_coverage": best["coverage_threshold"] if best else 0,
        "best_rank": best["rank"] if best else None,
        "top3_max_coverage": max((r["coverage_threshold"] for r in rows[:3]), default=0),
        "top5_max_coverage": max((r["coverage_threshold"] for r in rows[:5]), default=0),
        "pocket_count": len(rows), "pockets": rows
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crystal", required=True)
    parser.add_argument("--ligand", required=True, help="PDB residue name, e.g. MK1")
    parser.add_argument("--ligand-chain")
    parser.add_argument("--ligand-resseq")
    parser.add_argument("--pockets-pdb", required=True)
    parser.add_argument("--pockets-json", required=True)
    parser.add_argument("--threshold", type=float, default=3.0)
    parser.add_argument("--output")
    args = parser.parse_args()
    if not math.isfinite(args.threshold) or args.threshold <= 0:
        parser.error("--threshold must be positive")
    result = evaluate(args.crystal, args.ligand, args.pockets_pdb, args.pockets_json,
                      args.ligand_chain, args.ligand_resseq, args.threshold)
    encoded = json.dumps(result, indent=2)
    if args.output:
        Path(args.output).write_text(encoded + "\n", encoding="utf-8")
    else:
        print(encoded)


if __name__ == "__main__":
    main()
