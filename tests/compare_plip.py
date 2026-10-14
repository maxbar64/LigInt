#!/usr/bin/env python3
"""Confronta ligint con PLIP su uno o piu' PDB (precisione/sensibilita' per tipo di interazione).

Uso:
    python3 tests/compare_plip.py 1AZM.pdb 3PTB.pdb ...            ligint con le sue soglie
    python3 tests/compare_plip.py --plip-thr *.pdb                   ligint con soglie "alla PLIP"
    python3 tests/compare_plip.py -v *.pdb                           elenca anche i singoli disaccordi
    python3 tests/compare_plip.py --first-copy *.pdb                 una sola copia per ligando e PDB
    python3 tests/compare_plip.py --sdf-dir DIR *.pdb                usa DIR/<pdb>_<LIG>.sdf se presente

Requisiti: pip install plip openbabel-wheel

Come si confronta: PLIP analizza da solo tutti i ligandi (piccole molecole, anche con ioni) del PDB; per ogni
sito si lancia ligint sullo stesso ligando. Le interazioni sono confrontate a livello di RESIDUO e TIPO
(hydrophobic, hbond, salt, water, pistack, picat, halogen): (tipo, residuo, numero, catena).

Attenzione: PLIP NON e' la verita'. Ha le sue convenzioni (protonazione con Open Babel, soglie diverse,
centroidi dei gruppi carichi). Un disaccordo e' un motivo per guardare il caso, non un errore automatico.
Non sono confrontati: clash, repulsioni, legami covalenti, contatti con la catena principale e coordinazione
metallica (PLIP la descrive tra proteina e metallo, ligint tra ligando e metallo).
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
LIG_TYPES = {"Contatto idrofobico": "hydrophobic", "Legame a idrogeno": "hbond", "Ponte d'acqua": "water",
             "Ponte salino": "salt", "pi-stacking (parallelo)": "pistack", "pi-stacking (a T)": "pistack",
             "Cation-pi": "picat", "Legame ad alogeno": "halogen"}
CATS = ["hydrophobic", "hbond", "water", "salt", "pistack", "picat", "halogen"]
# solo le soglie confrontabili: il limite di 4,1 A degli H-bond di PLIP vale con la geometria degli H, ligint senza H no
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
                raise RuntimeError("PLIP non ha prodotto il report per %s:\n%s" % (pdb, p.stderr[-500:]))
            xml = os.path.join(outdir, cand[0])
    return xml


def parse_plip(xml):
    sites = []
    for bs in ET.parse(xml).getroot().findall("bindingsite"):
        ident = bs.find("identifiers")
        ligtype = ident.findtext("ligtype", "")
        if not ligtype.startswith("SMALLMOLECULE"):
            continue                                    # peptidi, DNA/RNA, ioni isolati: fuori da questo confronto
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
        raise RuntimeError("ligint fallito su %s %s: %s" % (pdb, site["hetid"], p.stderr[-300:]))
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
        if cat == "hbond" and "assistito da carica" in f[6]:   # H-bond tra gruppi di carica opposta: come il ponte salino, fuori dal punteggio
            continue
        inter.add((cat, res, num, chain))
        detail.setdefault((cat, res, num, chain), []).append("%s -> %s d=%s %s" % (f[3].split(":")[-1], f[4].split(":")[-1], f[5], f[6]))
        # ligint fonde ponte salino e H-bond sullo stesso contatto; PLIP li riporta insieme solo nel ~13% dei casi:
        # questi H-bond "derivati" NON entrano nel punteggio dell'H-bond (convenzione di PLIP)
    return inter, detail


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", nargs="+")
    ap.add_argument("--ligint", default="", help="binario di ligint (default: compila i sorgenti)")
    ap.add_argument("--plip-thr", action="store_true", help="soglie di ligint simili a quelle di PLIP")
    ap.add_argument("--plip-dir", default="", help="cartella per i report PLIP (riusati se gia' presenti)")
    ap.add_argument("--sdf-dir", default="", help="cartella con <pdb>_<LIG>.sdf (facoltativo)")
    ap.add_argument("--extra", default="", help="opzioni aggiuntive per ligint, tra apici")
    ap.add_argument("--first-copy", action="store_true",
                    help="conta una sola copia per (PDB, ligando): evita che i complessi con molte catene dominino i totali")
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
    print("%-28s %5s | %s" % ("sito", "PLIP", "  ".join("%-11s" % c for c in CATS)))
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
                        print("      solo ligint: %-11s %s %s%s  [%s]" % (c, x[1], x[2], x[3], "; ".join(ldet.get(x, []))[:90]))
                    for x in sorted(p - l):
                        print("      solo PLIP:   %-11s %s %s%s  [%s]" % (c, x[1], x[2], x[3], "; ".join(s["detail"].get(x, []))[:90]))
            print("%-28s %5d | %s" % ("%s %s:%s:%s" % (base, s["hetid"], s["chain"], s["pos"]), len(s["inter"]), "  ".join(cells)))
    print("\n%-12s %5s %5s %5s  %12s %12s" % ("tipo", "TP", "FP", "FN", "sensibilita'", "precisione"))
    gtp = gfp = gfn = 0
    for c in CATS:
        t = tot[c]
        if not (t["tp"] or t["fp"] or t["fn"]):
            continue
        rec = t["tp"] / (t["tp"] + t["fn"]) if (t["tp"] + t["fn"]) else float("nan")
        pre = t["tp"] / (t["tp"] + t["fp"]) if (t["tp"] + t["fp"]) else float("nan")
        print("%-12s %5d %5d %5d  %11.0f%% %11.0f%%" % (c, t["tp"], t["fp"], t["fn"], 100 * rec, 100 * pre))
        gtp += t["tp"]; gfp += t["fp"]; gfn += t["fn"]
    print("%-12s %5d %5d %5d  %11.0f%% %11.0f%%" % ("TOTALE", gtp, gfp, gfn, 100 * gtp / max(1, gtp + gfn),
                                                    100 * gtp / max(1, gtp + gfp)))
    print("(sensibilita' = quante interazioni PLIP ritrova ligint; precisione = quante di ligint sono anche in PLIP)")


if __name__ == "__main__":
    main()
