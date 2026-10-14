#!/usr/bin/env python3
"""Automatic benchmark of ligint against PLIP: from the PDB files to the report, in a single command.

Usage:
    python3 tests/benchmark.py archive.tar.gz                      an archive with the PDB files
    python3 tests/benchmark.py pdb_folder/                         a folder (also several arguments, also single files)
    python3 tests/benchmark.py pdb/ --groups groups.txt            separate report per group (e.g. A = tuning, B = blind)
    python3 tests/benchmark.py pdb/ --out results/ --all-copies    all the copies of the ligands (default: one per ligand)
    python3 tests/benchmark.py pdb/ "--extra=--plip-like"          options for ligint (values starting with "--" need the "=")

Groups file: one line "ID group" (e.g. "1ATP A"); IDs not listed go in the group "others".

Writes in --out (default ./benchmark_out):
    report.md           sensitivity/precision tables per type and group, causes of the disagreements, outcome per structure
    disagreements.csv   every disagreement with atoms, distances and probable cause (to filter/sort)
    per_category.csv    TP/FP/FN per site and type (for bootstrap and paired comparisons between configurations)
    plip/               PLIP XML reports (reused in later runs: PDB files already done are not processed again)

Requirements: pip install plip openbabel-wheel   (a C++17 compiler to build ligint, or --ligint BINARY)

How to read the result: PLIP is not the truth. The "causes" are automatic classifications of convenience
(distance beyond the ligint threshold, histidine partner, etc.): they indicate where to look, not who is right.
"""
import argparse
import collections
import csv
import datetime
import os
import re
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_plip as cp   # noqa: E402

MODES = [("ligint thresholds", []), ("PLIP-like thresholds (salt bridges 5.5 A, halogens 4.0 A)", cp.PLIP_LIKE)]


def gather(inputs, tmp):
    pdbs = []
    for p in inputs:
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                pdbs += [os.path.join(root, f) for f in sorted(files) if f.lower().endswith((".pdb", ".ent"))]
        elif os.path.isfile(p) and tarfile.is_tarfile(p):      # archive recognized by its content, not by its name
            d = tempfile.mkdtemp(prefix="pdbs_", dir=tmp)
            with tarfile.open(p) as t:
                t.extractall(d)
            for root, _, files in os.walk(d):
                pdbs += [os.path.join(root, f) for f in sorted(files) if f.lower().endswith((".pdb", ".ent"))]
        elif os.path.isfile(p):
            pdbs.append(p)
        else:
            sys.exit("not found: " + p)
    seen, out = set(), []
    for f in pdbs:
        k = os.path.splitext(os.path.basename(f))[0].upper()
        if k not in seen:
            seen.add(k)
            out.append((k, f))
    return out


def dist_of(detail):
    m = re.search(r"d=([0-9.]+)", detail or "")
    return float(m.group(1)) if m else None


def mechanism(cat, who, res, detail):
    """Probable cause of a disagreement (automatic classification of convenience)."""
    d = dist_of(detail)
    if cat == "hbond":
        if who == "PLIP":
            if d is None:
                return "PLIP: distance not available"
            return "PLIP: D-A beyond 3.5 A (outside the ligint threshold)" if d > 3.5 else "PLIP: D-A within 3.5 A, to be examined"
        if res == "HIS":
            return "ligint: histidine partner (tautomer assigned by PLIP)"
        if d is not None and d <= 3.0:
            return "ligint: strong contact (<=3.0 A) not reported by PLIP"
        return "ligint: 3.0-3.5 A contact not reported by PLIP"
    if cat == "salt":
        if who == "PLIP":
            return "PLIP: beyond 4.0 A (outside the ligint threshold)" if (d or 0) > 4.0 else "PLIP: within 4.0 A, to be examined"
        return "ligint: salt bridge not reported by PLIP"
    if cat == "water":
        return "ligint: water bridge not reported by PLIP (PLIP applies its own filters)" if who == "ligint" \
            else "PLIP: water bridge not reported by ligint"
    if cat == "halogen":
        return "ligint: halogen bond not reported by PLIP" if who == "ligint" else "PLIP: halogen bond not reported by ligint (F not handled, threshold)"
    if cat == "picat":
        return "ligint: cation-pi not reported by PLIP (ligand protonation?)" if who == "ligint" else "PLIP: cation-pi not reported by ligint"
    if cat == "pistack":
        return "ligint: stacking not reported by PLIP" if who == "ligint" else "PLIP: stacking not reported by ligint"
    return "ligint: hydrophobic contact not reported by PLIP" if who == "ligint" else "PLIP: hydrophobic contact not reported by ligint"


def pct(a, b):
    return "%3.0f%%" % (100.0 * a / b) if b else "  - "


def table(counts):
    rows = ["| type | TP | FP | FN | sensitivity | precision |", "|---|---:|---:|---:|---:|---:|"]
    gtp = gfp = gfn = 0
    for c in cp.CATS:
        t = counts[c]
        if not (t["tp"] or t["fp"] or t["fn"]):
            continue
        rows.append("| %s | %d | %d | %d | %s | %s |" % (c, t["tp"], t["fp"], t["fn"], pct(t["tp"], t["tp"] + t["fn"]), pct(t["tp"], t["tp"] + t["fp"])))
        gtp += t["tp"]; gfp += t["fp"]; gfn += t["fn"]
    rows.append("| **total** | %d | %d | %d | %s | %s |" % (gtp, gfp, gfn, pct(gtp, gtp + gfn), pct(gtp, gtp + gfp)))
    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="archive (.tar.gz), folder or PDB file")
    ap.add_argument("--out", default="benchmark_out")
    ap.add_argument("--groups", default="", help="file 'ID group' for the separate report")
    ap.add_argument("--all-copies", action="store_true", help="count all the copies of each ligand")
    ap.add_argument("--ligint", default="", help="ligint binary (default: build from the sources)")
    ap.add_argument("--extra", default="", help="additional options for ligint, in quotes")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="bench_")
    binary = a.ligint
    if not binary:
        binary = os.path.join(tmp, "ligint")
        cp_build = __import__("_build")
        cp_build.build_ligint(binary)
    plipdir = os.path.join(a.out, "plip")
    os.makedirs(plipdir, exist_ok=True)
    groups = {}
    if a.groups:
        for line in open(a.groups):
            f = line.split()
            if len(f) >= 2:
                groups[f[0].upper()] = f[1]
    pdbs = gather(a.inputs, tmp)
    print("%d PDB files to process" % len(pdbs))

    counts = {(m, g): {c: collections.Counter() for c in cp.CATS} for m, _ in MODES for g in set(groups.get(k, "others") for k, _ in pdbs) | {"ALL"}}
    disagreements, per_struct, per_cat, failures = [], [], [], []
    for pid, path in pdbs:
        grp = groups.get(pid, "others")
        try:
            sites = cp.parse_plip(cp.run_plip(path, plipdir))
        except Exception as e:                       # noqa: BLE001
            failures.append((pid, "PLIP: %s" % str(e)[:200])); print("  %s: PLIP failed" % pid); continue
        if not sites:
            failures.append((pid, "no small-molecule ligand for PLIP")); continue
        seen = set()
        for s in sites:
            if not a.all_copies:
                if s["hetid"] in seen:
                    continue
                seen.add(s["hetid"])
            for mi, (mname, mopts) in enumerate(MODES):
                try:
                    li, ldet = cp.run_ligint(binary, path, s, mopts + a.extra.split())
                except Exception as e:               # noqa: BLE001
                    failures.append((pid, "ligint on %s: %s" % (s["hetid"], str(e)[:160]))); continue
                tot = collections.Counter()
                for c in cp.CATS:
                    p = {x for x in s["inter"] if x[0] == c}
                    l = {x for x in li if x[0] == c}
                    tp, fp, fn = len(p & l), len(l - p), len(p - l)
                    for g in (grp, "ALL"):
                        counts[(mname, g)][c].update(tp=tp, fp=fp, fn=fn)
                    tot.update(tp=tp, fp=fp, fn=fn)
                    if mi == 0:
                        per_cat.append([pid, grp, "%s:%s:%s" % (s["hetid"], s["chain"], s["pos"]), c, tp, fp, fn])
                    if mi == 0:                      # list of the disagreements with the ligint thresholds
                        for x in sorted(l - p):
                            det = "; ".join(ldet.get(x, []))
                            disagreements.append([pid, grp, "%s:%s:%s" % (s["hetid"], s["chain"], s["pos"]), c, "ligint",
                                                  "%s %s%s" % (x[1], x[2], x[3]), det, mechanism(c, "ligint", x[1], det)])
                        for x in sorted(p - l):
                            det = "; ".join(s["detail"].get(x, []))
                            disagreements.append([pid, grp, "%s:%s:%s" % (s["hetid"], s["chain"], s["pos"]), c, "PLIP",
                                                  "%s %s%s" % (x[1], x[2], x[3]), det, mechanism(c, "PLIP", x[1], det)])
                if mi == 0:
                    per_struct.append((pid, grp, "%s:%s:%s" % (s["hetid"], s["chain"], s["pos"]), len(s["inter"]),
                                       tot["tp"], tot["fp"], tot["fn"]))
        print("  %s ok" % pid)

    with open(os.path.join(a.out, "disagreements.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdb", "group", "site", "type", "only_in", "residue", "detail", "probable_cause"])
        w.writerows(disagreements)

    with open(os.path.join(a.out, "per_category.csv"), "w", newline="") as f:      # counts per site and type (for bootstrap and paired comparisons)
        w = csv.writer(f)
        w.writerow(["pdb", "group", "site", "type", "tp", "fp", "fn"])
        w.writerows(per_cat)

    out = []
    out.append("# Benchmark of ligint against PLIP\n")
    out.append("Date: %s. PDB files: %d, sites compared: %d%s.\n" % (datetime.date.today().isoformat(), len(pdbs), len(per_struct),
                                                                  " (all copies)" if a.all_copies else " (one copy per ligand and PDB)"))
    out.append("Sensitivity = how many PLIP interactions ligint finds. Precision = how many of ligint's are also in PLIP. "
               "PLIP is not the truth: a disagreement is a reason to look at the case.\n")
    gnames = sorted({g for (_, g) in counts if g != "ALL"}) + ["ALL"]
    for mname, _ in MODES:
        out.append("## %s\n" % mname)
        for g in gnames:
            if any(counts[(mname, g)][c] for c in cp.CATS):
                out.append("### group %s\n" % g)
                out.append(table(counts[(mname, g)]) + "\n")
    mech = collections.Counter((d[3], d[7]) for d in disagreements)
    out.append("## Probable causes of the disagreements (ligint thresholds, all groups)\n")
    out.append("| type | cause | n |\n|---|---|---:|")
    for (c, m), n in sorted(mech.items(), key=lambda kv: (cp.CATS.index(kv[0][0]), -kv[1])):
        out.append("| %s | %s | %d |" % (c, m, n))
    out.append("")
    out.append("## Outcome per structure (ligint thresholds)\n")
    out.append("| PDB | group | site | PLIP interactions | TP | FP | FN |\n|---|---|---|---:|---:|---:|---:|")
    for r in per_struct:
        out.append("| %s | %s | %s | %d | %d | %d | %d |" % r)
    if failures:
        out.append("\n## Problems\n")
        for pid, msg in failures:
            out.append("- %s: %s" % (pid, msg))
    with open(os.path.join(a.out, "report.md"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("\nwritten: %s/report.md and %s/disagreements.csv" % (a.out, a.out))
    for g in gnames:
        c = counts[(MODES[0][0], g)]
        tp = sum(c[k]["tp"] for k in cp.CATS); fp = sum(c[k]["fp"] for k in cp.CATS); fn = sum(c[k]["fn"] for k in cp.CATS)
        if tp + fp + fn:
            print("group %-6s sensitivity %s  precision %s" % (g, pct(tp, tp + fn).strip(), pct(tp, tp + fp).strip()))


if __name__ == "__main__":
    main()
