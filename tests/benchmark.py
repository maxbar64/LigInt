#!/usr/bin/env python3
"""Benchmark automatico di ligint contro PLIP: dai PDB al report, in un solo comando.

Uso:
    python3 tests/benchmark.py archivio.tar.gz                     un archivio con i PDB
    python3 tests/benchmark.py cartella_pdb/                       una cartella (anche piu' argomenti, anche singoli file)
    python3 tests/benchmark.py pdb/ --groups gruppi.txt            report separato per gruppo (es. A = taratura, B = cieco)
    python3 tests/benchmark.py pdb/ --out risultati/ --all-copies  tutte le copie dei ligandi (default: una per ligando)
    python3 tests/benchmark.py pdb/ "--extra=--plip-like"           opzioni per ligint (con valori che iniziano per "--" serve il "=")

File dei gruppi: una riga "ID gruppo" (es. "1ATP A"); gli ID non elencati vanno nel gruppo "altri".

Produce in --out (default ./benchmark_out):
    report.md           tabelle di sensibilita'/precisione per tipo e gruppo, cause dei disaccordi, esito per struttura
    disagreements.csv   ogni disaccordo con atomi, distanze e causa probabile (da filtrare/ordinare)
    per_category.csv    TP/FP/FN per sito e tipo (per bootstrap e confronti appaiati tra configurazioni)
    plip/               report XML di PLIP (riusati alle esecuzioni successive: i PDB gia' fatti non si rielaborano)

Requisiti: pip install plip openbabel-wheel   (un compilatore C++17 per costruire ligint, oppure --ligint BINARIO)

Come leggere il risultato: PLIP non e' la verita'. Le "cause" sono classificazioni automatiche di comodo
(distanza oltre la soglia di ligint, partner istidina, ecc.): indicano dove guardare, non chi ha ragione.
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

MODES = [("soglie di ligint", []), ("soglie alla PLIP (ponti salini 5,5 A, alogeni 4,0 A)", cp.PLIP_LIKE)]


def gather(inputs, tmp):
    pdbs = []
    for p in inputs:
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                pdbs += [os.path.join(root, f) for f in sorted(files) if f.lower().endswith((".pdb", ".ent"))]
        elif os.path.isfile(p) and tarfile.is_tarfile(p):      # archivio riconosciuto dal contenuto, non dal nome
            d = tempfile.mkdtemp(prefix="pdbs_", dir=tmp)
            with tarfile.open(p) as t:
                t.extractall(d)
            for root, _, files in os.walk(d):
                pdbs += [os.path.join(root, f) for f in sorted(files) if f.lower().endswith((".pdb", ".ent"))]
        elif os.path.isfile(p):
            pdbs.append(p)
        else:
            sys.exit("non trovato: " + p)
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
    """Causa probabile di un disaccordo (classificazione automatica di comodo)."""
    d = dist_of(detail)
    if cat == "hbond":
        if who == "PLIP":
            if d is None:
                return "PLIP: distanza non disponibile"
            return "PLIP: D-A oltre 3,5 A (fuori soglia di ligint)" if d > 3.5 else "PLIP: D-A entro 3,5 A, da esaminare"
        if res == "HIS":
            return "ligint: partner istidina (tautomero assegnato da PLIP)"
        if d is not None and d <= 3.0:
            return "ligint: contatto forte (<=3,0 A) non riportato da PLIP"
        return "ligint: contatto 3,0-3,5 A non riportato da PLIP"
    if cat == "salt":
        if who == "PLIP":
            return "PLIP: oltre 4,0 A (fuori soglia di ligint)" if (d or 0) > 4.0 else "PLIP: entro 4,0 A, da esaminare"
        return "ligint: ponte salino non riportato da PLIP"
    if cat == "water":
        return "ligint: ponte d'acqua non riportato da PLIP (filtri propri di PLIP)" if who == "ligint" \
            else "PLIP: ponte d'acqua non riportato da ligint"
    if cat == "halogen":
        return "ligint: alogeno non riportato" if who == "ligint" else "PLIP: alogeno non riportato da ligint (F non trattato, soglia)"
    if cat == "picat":
        return "ligint: cation-pi non riportato da PLIP (protonazione del ligando?)" if who == "ligint" else "PLIP: cation-pi non riportato"
    if cat == "pistack":
        return "ligint: stacking non riportato da PLIP" if who == "ligint" else "PLIP: stacking non riportato da ligint"
    return "ligint: contatto idrofobico non riportato da PLIP" if who == "ligint" else "PLIP: contatto idrofobico non riportato da ligint"


def pct(a, b):
    return "%3.0f%%" % (100.0 * a / b) if b else "  - "


def table(counts):
    rows = ["| tipo | TP | FP | FN | sensibilità | precisione |", "|---|---:|---:|---:|---:|---:|"]
    gtp = gfp = gfn = 0
    for c in cp.CATS:
        t = counts[c]
        if not (t["tp"] or t["fp"] or t["fn"]):
            continue
        rows.append("| %s | %d | %d | %d | %s | %s |" % (c, t["tp"], t["fp"], t["fn"], pct(t["tp"], t["tp"] + t["fn"]), pct(t["tp"], t["tp"] + t["fp"])))
        gtp += t["tp"]; gfp += t["fp"]; gfn += t["fn"]
    rows.append("| **totale** | %d | %d | %d | %s | %s |" % (gtp, gfp, gfn, pct(gtp, gtp + gfn), pct(gtp, gtp + gfp)))
    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="archivio .tar.gz, cartella o file PDB")
    ap.add_argument("--out", default="benchmark_out")
    ap.add_argument("--groups", default="", help="file 'ID gruppo' per il report separato")
    ap.add_argument("--all-copies", action="store_true", help="conta tutte le copie di ogni ligando")
    ap.add_argument("--ligint", default="", help="binario di ligint (default: compila i sorgenti)")
    ap.add_argument("--extra", default="", help="opzioni aggiuntive per ligint, tra apici")
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
    print("%d PDB da elaborare" % len(pdbs))

    counts = {(m, g): {c: collections.Counter() for c in cp.CATS} for m, _ in MODES for g in set(groups.get(k, "altri") for k, _ in pdbs) | {"TUTTI"}}
    disagreements, per_struct, per_cat, failures = [], [], [], []
    for pid, path in pdbs:
        grp = groups.get(pid, "altri")
        try:
            sites = cp.parse_plip(cp.run_plip(path, plipdir))
        except Exception as e:                       # noqa: BLE001
            failures.append((pid, "PLIP: %s" % str(e)[:200])); print("  %s: PLIP fallito" % pid); continue
        if not sites:
            failures.append((pid, "nessun ligando a piccola molecola per PLIP")); continue
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
                    failures.append((pid, "ligint su %s: %s" % (s["hetid"], str(e)[:160]))); continue
                tot = collections.Counter()
                for c in cp.CATS:
                    p = {x for x in s["inter"] if x[0] == c}
                    l = {x for x in li if x[0] == c}
                    tp, fp, fn = len(p & l), len(l - p), len(p - l)
                    for g in (grp, "TUTTI"):
                        counts[(mname, g)][c].update(tp=tp, fp=fp, fn=fn)
                    tot.update(tp=tp, fp=fp, fn=fn)
                    if mi == 0:
                        per_cat.append([pid, grp, "%s:%s:%s" % (s["hetid"], s["chain"], s["pos"]), c, tp, fp, fn])
                    if mi == 0:                      # elenco dei disaccordi con le soglie di ligint
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
        w.writerow(["pdb", "gruppo", "sito", "tipo", "solo_in", "residuo", "dettaglio", "causa_probabile"])
        w.writerows(disagreements)

    with open(os.path.join(a.out, "per_category.csv"), "w", newline="") as f:      # conteggi per sito e tipo (per bootstrap e confronti appaiati)
        w = csv.writer(f)
        w.writerow(["pdb", "gruppo", "sito", "tipo", "tp", "fp", "fn"])
        w.writerows(per_cat)

    out = []
    out.append("# Benchmark ligint contro PLIP\n")
    out.append("Data: %s. PDB: %d, siti confrontati: %d%s.\n" % (datetime.date.today().isoformat(), len(pdbs), len(per_struct),
                                                                  " (tutte le copie)" if a.all_copies else " (una copia per ligando e PDB)"))
    out.append("Sensibilità = quante interazioni di PLIP ritrova ligint. Precisione = quante di ligint sono anche in PLIP. "
               "PLIP non è la verità: un disaccordo è un motivo per guardare il caso.\n")
    gnames = sorted({g for (_, g) in counts if g != "TUTTI"}) + ["TUTTI"]
    for mname, _ in MODES:
        out.append("## %s\n" % mname)
        for g in gnames:
            if any(counts[(mname, g)][c] for c in cp.CATS):
                out.append("### gruppo %s\n" % g)
                out.append(table(counts[(mname, g)]) + "\n")
    mech = collections.Counter((d[3], d[7]) for d in disagreements)
    out.append("## Cause probabili dei disaccordi (soglie di ligint, tutti i gruppi)\n")
    out.append("| tipo | causa | n |\n|---|---|---:|")
    for (c, m), n in sorted(mech.items(), key=lambda kv: (cp.CATS.index(kv[0][0]), -kv[1])):
        out.append("| %s | %s | %d |" % (c, m, n))
    out.append("")
    out.append("## Esito per struttura (soglie di ligint)\n")
    out.append("| PDB | gruppo | sito | interazioni PLIP | TP | FP | FN |\n|---|---|---|---:|---:|---:|---:|")
    for r in per_struct:
        out.append("| %s | %s | %s | %d | %d | %d | %d |" % r)
    if failures:
        out.append("\n## Problemi\n")
        for pid, msg in failures:
            out.append("- %s: %s" % (pid, msg))
    with open(os.path.join(a.out, "report.md"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("\nscritti: %s/report.md e %s/disagreements.csv" % (a.out, a.out))
    for g in gnames:
        c = counts[(MODES[0][0], g)]
        tp = sum(c[k]["tp"] for k in cp.CATS); fp = sum(c[k]["fp"] for k in cp.CATS); fn = sum(c[k]["fn"] for k in cp.CATS)
        if tp + fp + fn:
            print("gruppo %-6s sensibilita' %s  precisione %s" % (g, pct(tp, tp + fn).strip(), pct(tp, tp + fp).strip()))


if __name__ == "__main__":
    main()
