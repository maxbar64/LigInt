#!/usr/bin/env python3
"""Reproducible external pocket benchmark; crystallographic ligands are evaluation-only."""
import argparse
import json
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

CASES={"1STP":"BTN","1IEP":"STI"}

def ligand_instances(pdb, residue):
    result=defaultdict(list)
    for line in pdb.read_text(encoding="utf-8",errors="replace").splitlines():
        if not line.startswith("HETATM") or line[17:20].strip()!=residue:
            continue
        alt=line[16:17]
        if alt not in (" ","A"):continue
        element=line[76:78].strip().upper()
        if element in ("H","D"):continue
        chain=line[21:22].strip()
        seq=line[22:26].strip()
        result[(chain,seq)].append(line)
    return {k:len(v) for k,v in sorted(result.items())}

def run(args):
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    summaries=[]
    for pdb_id,ligand in CASES.items():
        pdb=out/(pdb_id+".pdb")
        if not pdb.exists():
            url=f"https://files.rcsb.org/download/{pdb_id}.pdb"
            urllib.request.urlretrieve(url,pdb)
        instances=ligand_instances(pdb,ligand)
        if not instances:raise RuntimeError(f"No {ligand} present in {pdb}")
        protein_chains={line[21:22].strip() for line in
            pdb.read_text(encoding="utf-8",errors="replace").splitlines()
            if line.startswith("ATOM  ")}
        for (chain,resseq),n_atoms in instances.items():
            if chain not in protein_chains:
                raise RuntimeError(f"Cannot match ligand {pdb_id} chain {chain!r} to protein chains")
            tag=f"{pdb_id}_{chain or 'blank'}_{resseq}"
            directory=out/tag
            directory.mkdir(exist_ok=True)
            # Select only protein chain corresponding to the ligand instance.
            # Other ligand instances do not contribute to pocket prediction.
            base=[str(args.binary),str(pdb),"--chains",chain,
                  "--rank-by","heuristic","--morphology","opening",
                  "--morph-radius","1",
                  "--output",str(directory/"opening.json"),
                  "--pdb",str(directory/"opening.pdb"),
                  "--groups",str(directory/"raw_groups.json"),
                  "--diverse-groups",str(directory/"diverse_groups.json")]
            subprocess.run(base,check=True)
            trials={}
            for mode in ("raw","diverse"):
                output=directory/(mode+"_benchmark.json")
                subprocess.run([sys.executable,str(args.benchmark_script),
                  "--crystal",str(pdb),"--ligand",ligand,
                  "--ligand-chain",chain,"--ligand-resseq",resseq,
                  "--pockets-pdb",str(directory/"opening.pdb"),
                  "--groups-json",str(directory/(mode+"_groups.json")),
                  "--distance","3.0","--output",str(output)],check=True)
                with output.open(encoding="utf-8") as f:
                    report=json.load(f)
                trials[mode]={"candidate_count":report["candidate_count"],
                    "first_full_coverage_rank":report["first_full_coverage_rank"],
                    "best_coverage":report["best_coverage"],
                    "top1_coverage":report["results"][0]["coverage"] if report["results"] else 0,
                    "top3_best_coverage":max((r["coverage"] for r in report["results"][:3]),default=0)}
            summaries.append({"pdb":pdb_id,"ligand":ligand,"ligand_chain":chain,
                "ligand_resseq":resseq,"ligand_heavy_atoms":n_atoms,
                "protein_chains_used":[chain],"raw":trials["raw"],"diverse":trials["diverse"]})
    with (out/"external_summary.json").open("w",encoding="utf-8") as f:
        json.dump({"settings":{"morphology":"opening","morph_radius":1,
            "rank_by":"heuristic","distance_a":3.0},
            "cases":summaries},f,indent=2)
        f.write("\n")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary",type=Path,required=True)
    p.add_argument("--benchmark-script",type=Path,default=Path(__file__).with_name("pocket_group_benchmark.py"))
    p.add_argument("--output-dir",type=Path,default=Path("pocket_external_results"))
    args=p.parse_args()
    args.binary=args.binary.resolve()
    args.benchmark_script=args.benchmark_script.resolve()
    args.output_dir=args.output_dir.resolve()
    run(args)
if __name__=="__main__":
    main()
