#!/usr/bin/env python3
"""Compare ligint with PLIP on one or more PDB files (precision/sensitivity per interaction type).

Usage:
    python3 tests/compare_plip.py 1AZM.pdb 3PTB.pdb ...            ligint with its own thresholds
    python3 tests/compare_plip.py --plip-thr *.pdb                   ligint with "PLIP-like" thresholds
    python3 tests/compare_plip.py -v *.pdb                           also list the individual disagreements
    python3 tests/compare_plip.py --first-copy *.pdb                 one copy per ligand and PDB
    python3 tests/compare_plip.py --sdf-dir DIR *.pdb                use DIR/<pdb>_<LIG>.sdf if present

Requirements: pip install plip openbabel-wheel

How the comparison works: PLIP analyzes all the ligands (small molecules, also with ions) of the PDB by itself; for
each site ligint is run on the same ligand. The interactions are compared at the level of RESIDUE and TYPE
(hydrophobic, hbond, salt, water, pistack, picat, halogen): (type, residue, number, chain).

Warning: PLIP is NOT the truth. It has its own conventions (protonation with Open Babel, different thresholds,
centroids of charged groups). A disagreement is a reason to look at the case, not an automatic error.
Not compared: clashes, repulsions, covalent bonds, contacts with the backbone and metal coordination
(PLIP describes it between protein and metal, ligint between ligand and metal).
"""
import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _build  # noqa: E402
PLIP_TAGS = {"hydrophobic_interaction": "hydrophobic", "hydrogen_bond": "hbond", "water_bridge": "water",
             "salt_bridge": "salt", "pi_stack": "pistack", "pi_cation_interaction": "picat",
             "halogen_bond": "halogen"}
LIG_TYPES = {"Hydrophobic contact": "hydrophobic", "Hydrogen bond": "hbond", "Water bridge": "water",
             "Salt bridge": "salt", "pi-stacking (parallel)": "pistack", "pi-stacking (T-shaped)": "pistack",
             "Cation-pi": "picat", "Halogen bond": "halogen"}
CATS = ["hydrophobic", "hbond", "water", "salt", "pistack", "picat", "halogen"]
# only the comparable thresholds: the 4.1 A limit of PLIP's H-bonds applies with computed H geometry, ligint without H does not
PLIP_LIKE = ["--salt-max", "5.5", "--hal-max", "4.0"]


def run_plip(pdb, outdir):
    xml = os.path.join(outdir, os.path.splitext(os.path.basename(pdb))[0].upper() + "_report.xml")
    if not os.path.exists(xml):
        p = subprocess.run([sys.executable, "-m", "plip.plipcmd", "-f", pdb, "-x", "-o", outdir],
                           capture_output=True, text=True)
        if not os.path.exists(xml):
            cand = [f for f in os.listdir(outdir) if f.endswith("_report.xml")
                    and f.upper().startswith(os.path.splitext(os.path.basename(pdb))[0].upper())]
            if not cand:
                raise RuntimeError("PLIP did not produce the report for %s:\n%s" % (pdb, p.stderr[-500:]))
            xml = os.path.join(outdir, cand[0])
    return xml


def parse_plip(xml):
    sites = []
    for bs in ET.parse(xml).getroot().findall("bindingsite"):
        ident = bs.find("identifiers")
        ligtype = ident.findtext("ligtype", "")
        if not ligtype.startswith("SMALLMOLECULE"):
            continue                                    # peptides, DNA/RNA, isolated ions: outside this comparison
        inter = set()
        detail = {}
        node = bs.find("interactions")
        if node is not None:
            for grp in node:
                for it in grp:
                    cat = PLIP_TAGS.get(it.tag)
                    if cat:
                        key = (cat, it.findtext("restype"), it.findtext("resnr"), it.findtext("reschain"))
                        inter.add(key)
                        d = it.findtext("dist") or it.findtext("dist_d-a") or it.findtext("distance") or it.findtext("centdist") or "?"
                        extra = {"hbond": " don=%s acc=%s" % (it.findtext("donortype"), it.findtext("acceptortype")),
                                 "pistack": " ang=%s off=%s" % (it.findtext("angle"), it.findtext("offset")),
                                 "picat": " off=%s" % it.findtext("offset")}.get(cat, "")
                        detail.setdefault(key, []).append("d=%s%s" % (d, extra))
        sites.append({"hetid": ident.findtext("hetid"), "chain": ident.findtext("chain"),
                      "pos": ident.findtext("position"), "name": ident.findtext("longname"), "inter": inter,
                      "detail": detail})
    return sites


PA_RE = re.compile(r"^(\S+)\s+(-?\d+)(.*?):")


def run_ligint(binary, pdb, site, extra):
    cmd = [binary, pdb, site["hetid"], site["chain"], site["pos"], "--csv"] + extra
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError("ligint failed on %s %s: %s" % (pdb, site["hetid"], p.stderr[-300:]))
    inter = set()
    detail = {}
    for line in p.stdout.splitlines()[1:]:
        f = line.split(",", 6)
        if len(f) < 7 or f[1] not in LIG_TYPES:
            continue
        m = PA_RE.match(f[4])
        if not m:
            continue
        res, num, tail = m.group(1), m.group(2), m.group(3)
        chain = tail[-1] if tail else " "
        cat = LIG_TYPES[f[1]]
        if cat == "hbond" and "charge-assisted" in f[6]:   # H-bond between oppositely charged groups: like the salt bridge, outside the score
            continue
        inter.add((cat, res, num, chain))
        detail.setdefault((cat, res, num, chain), []).append("%s -> %s d=%s %s" % (f[3].split(":")[-1], f[4].split(":")[-1], f[5], f[6]))
        # ligint merges salt bridge and H-bond on the same contact; PLIP reports them together only in ~13% of cases:
        # these "derived" H-bonds do NOT enter the H-bond score (PLIP convention)
    return inter, detail


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", nargs="+")
    ap.add_argument("--ligint", default="", help="ligint binary (default: build from the sources)")
    ap.add_argument("--plip-thr", action="store_true", help="ligint thresholds similar to those of PLIP")
    ap.add_argument("--plip-dir", default="", help="folder for the PLIP reports (reused if already present)")
    ap.add_argument("--sdf-dir", default="", help="folder with <pdb>_<LIG>.sdf (optional)")
    ap.add_argument("--extra", default="", help="additional options for ligint, in quotes")
    ap.add_argument("--first-copy", action="store_true",
                    help="count only one copy per (PDB, ligand): prevents complexes with many chains from dominating the totals")
    ap.add_argument("-v", action="store_true")
    a = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="cmp_")
    binary = a.ligint
    if not binary:
        binary = os.path.join(tmp, "ligint")
        _build.build_ligint(binary)
    plipdir = a.plip_dir or os.path.join(tmp, "plip")
    os.makedirs(plipdir, exist_ok=True)
    extra = (PLIP_LIKE if a.plip_thr else []) + a.extra.split()

    tot = {c: collections.Counter() for c in CATS}
    print("%-28s %5s | %s" % ("site", "PLIP", "  ".join("%-11s" % c for c in CATS)))
    print("%-28s %5s | %s" % ("", "", "  ".join("%-11s" % "TP/FP/FN" for _ in CATS)))
    for pdb in a.pdb:
        base = os.path.splitext(os.path.basename(pdb))[0]
        sites = parse_plip(run_plip(pdb, plipdir))
        seen = set()
        for s in sites:
            if a.first_copy:
                if s["hetid"] in seen:
                    continue
                seen.add(s["hetid"])
            ex = list(extra)
            if a.sdf_dir:
                f = os.path.join(a.sdf_dir, "%s_%s.sdf" % (base, s["hetid"]))
                if os.path.exists(f):
                    ex += ["--sdf", f]
            li, ldet = run_ligint(binary, pdb, s, ex)
            cells = []
            for c in CATS:
                p = {x for x in s["inter"] if x[0] == c}
                l = {x for x in li if x[0] == c}
                tp, fp, fn = len(p & l), len(l - p), len(p - l)
                tot[c].update(tp=tp, fp=fp, fn=fn)
                cells.append("%-11s" % ("%d/%d/%d" % (tp, fp, fn) if (tp or fp or fn) else "."))
                if a.v and (fp or fn):
                    for x in sorted(l - p):
                        print("      ligint only: %-11s %s %s%s  [%s]" % (c, x[1], x[2], x[3], "; ".join(ldet.get(x, []))[:90]))
                    for x in sorted(p - l):
                        print("      PLIP only:   %-11s %s %s%s  [%s]" % (c, x[1], x[2], x[3], "; ".join(s["detail"].get(x, []))[:90]))
            print("%-28s %5d | %s" % ("%s %s:%s:%s" % (base, s["hetid"], s["chain"], s["pos"]), len(s["inter"]), "  ".join(cells)))
    print("\n%-12s %5s %5s %5s  %12s %12s" % ("type", "TP", "FP", "FN", "sensitivity", "precision"))
    gtp = gfp = gfn = 0
    for c in CATS:
        t = tot[c]
        if not (t["tp"] or t["fp"] or t["fn"]):
            continue
        rec = t["tp"] / (t["tp"] + t["fn"]) if (t["tp"] + t["fn"]) else float("nan")
        pre = t["tp"] / (t["tp"] + t["fp"]) if (t["tp"] + t["fp"]) else float("nan")
        print("%-12s %5d %5d %5d  %11.0f%% %11.0f%%" % (c, t["tp"], t["fp"], t["fn"], 100 * rec, 100 * pre))
        gtp += t["tp"]; gfp += t["fp"]; gfn += t["fn"]
    print("%-12s %5d %5d %5d  %11.0f%% %11.0f%%" % ("TOTAL", gtp, gfp, gfn, 100 * gtp / max(1, gtp + gfn),
                                                    100 * gtp / max(1, gtp + gfp)))
    print("(sensitivity = how many PLIP interactions ligint finds; precision = how many of ligint's are also in PLIP)")


if __name__ == "__main__":
    main()
