#!/usr/bin/env python3
"""GPCRdb reference -> PDB receptor-chain mapping. Standard library only.

Uses the GPCRdb REST residues/extended endpoint. Never invent generic positions.
Output is compatible with LigInt --gpcrdb-mapping (chain,resi,aa,gpcrdb,segment).
"""
import argparse
import json
import pathlib
import sys
from urllib.request import Request, urlopen

AA = dict(zip(
    "ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCEQGHILKMFPSTWYV"
))
AA.update({"MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C"})

def pdb_residues(path, chain):
    rows, seen = [], set()
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        if len(line) < 54 or not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[21].strip() != chain or line[12:16].strip() != "CA":
            continue
        if line[16] not in " A":
            continue
        code = line[17:20].strip().upper()
        aa = AA.get(code)
        if not aa:
            continue
        resi = line[22:26].strip()
        icode = line[26].strip()
        key = resi + icode
        if key in seen:
            continue
        seen.add(key)
        rows.append((resi, icode, aa, code))
    if not rows:
        raise ValueError(f"No receptor C-alpha residues found for chain {chain!r}")
    return rows

def identify_entry_from_pdb(pdb_id):
    """Best-effort GPCRdb structure lookup; fail closed for uncertain responses."""
    import re
    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        raise ValueError("PDB identifier must be four alphanumeric characters")
    url = f"https://gpcrdb.org/services/structure/{pdb_id.upper()}/"
    with urlopen(Request(url, headers={"User-Agent": "LigInt-GPCRdb/1.0",
                                        "Accept": "application/json"}), timeout=15) as res:
        data = json.load(res)
    if isinstance(data, list):
        if len(data) != 1:
            raise ValueError("GPCRdb structure lookup returned multiple or no records")
        data = data[0]
    if not isinstance(data, dict):
        raise ValueError("Unexpected GPCRdb structure response")
    candidates = [data.get("protein"), data.get("receptor"), data.get("protein_entry_name")]
    for value in candidates:
        if isinstance(value, dict):
            value = value.get("entry_name")
        if isinstance(value, str) and re.fullmatch(r"[a-z0-9_]+", value.lower()):
            return value.lower()
    raise ValueError("GPCRdb structure lookup does not contain a validated receptor entry")

def fetch_reference(entry, cache_dir):
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / (entry + ".residues.json")
    if target.exists():
        data = json.loads(target.read_text())
    else:
        url = f"https://gpcrdb.org/services/residues/extended/{entry}/"
        with urlopen(Request(url, headers={"User-Agent": "LigInt-GPCRdb/1.0", "Accept": "application/json"}), timeout=25) as res:
            data = json.load(res)
        if not isinstance(data, list) or not data:
            raise ValueError("Unexpected GPCRdb residues response")
        target.write_text(json.dumps(data, indent=2))
    return data

def gpcr_number(raw):
    if not raw:
        return ""
    raw = str(raw)
    if "x" in raw:
        prefix, after = raw.split("x", 1)
        return prefix.split(".")[0] + "x" + after
    return ""

def align(query, reference):
    """Semiglobal alignment: query chain must align; free reference end gaps.
    Pure Python affine-free Needleman-Wunsch; no gap/no mismatch is fabricated.
    """
    m, n = len(query), len(reference)
    if m*n > 5_000_000:
        raise ValueError("Sequence alignment exceeds safety limit")
    score = [[0]*(n+1) for _ in range(m+1)]
    trace = [[0]*(n+1) for _ in range(m+1)]
    for i in range(1, m+1):
        score[i][0] = -4*i
        trace[i][0] = 1
    for i in range(1, m+1):
        qi = query[i-1]
        for j in range(1, n+1):
            diag = score[i-1][j-1] + (3 if qi == reference[j-1] else -3)
            up = score[i-1][j] - 4
            left = score[i][j-1] - 4
            if diag >= up and diag >= left: score[i][j], trace[i][j] = diag, 0
            elif up >= left: score[i][j], trace[i][j] = up, 1
            else: score[i][j], trace[i][j] = left, 2
    j = max(range(n+1), key=lambda k: score[m][k])
    i, pairs = m, []
    while i:
        if j and trace[i][j] == 0:
            pairs.append((i-1, j-1))
            i -= 1
            j -= 1
        elif trace[i][j] == 1 or j == 0:
            i -= 1
        else:
            j -= 1
    return list(reversed(pairs))

def detect_chain(pdb, reference):
    """Identify the best receptor chain by alignment, rejecting ambiguous matches."""
    chains = sorted({
        line[21].strip() for line in pathlib.Path(pdb).read_text(errors="replace").splitlines()
        if len(line) >= 27 and line.startswith("ATOM  ") and line[12:16].strip() == "CA"
        and line[21].strip()
    })
    candidates = []
    target = "".join(r[1] for r in reference)
    for chain in chains:
        observed = pdb_residues(pdb, chain)
        if len(observed) < 35:
            continue
        pairs = align("".join(r[2] for r in observed), target)
        matches = sum(observed[i][2] == reference[j][1] for i, j in pairs)
        identity = matches / len(pairs) if pairs else 0
        coverage = matches / len(observed)
        if matches >= 30 and identity >= .85 and coverage >= .65:
            candidates.append((matches, identity, coverage, chain))
    candidates.sort(reverse=True)
    if not candidates:
        raise ValueError("No chain matches this GPCRdb receptor. Select the chain and check the receptor entry.")
    if len(candidates) > 1 and candidates[1][0] >= candidates[0][0] * .95:
        raise ValueError("Multiple similarly matching receptor chains: " +
                         ", ".join(item[3] for item in candidates[:3]) +
                         ". Select one explicitly.")
    return candidates[0][3]

def create(pdb, chain, entry, cache_dir, output, audit_path):
    raw = fetch_reference(entry, cache_dir)
    reference = []
    for item in raw:
        number = item.get("sequence_number")
        aa = str(item.get("amino_acid") or "").upper()
        if number is None or len(aa) != 1:
            continue
        reference.append((int(number), aa, gpcr_number(item.get("display_generic_number")),
                          str(item.get("protein_segment") or "")))
    reference.sort(key=lambda r: r[0])
    if len(reference) < 50:
        raise ValueError("Reference has too few GPCRdb residues")
    if chain.lower() == "auto":
        chain = detect_chain(pdb, reference)
    observed = pdb_residues(pdb, chain)
    pairs = align("".join(x[2] for x in observed), "".join(x[1] for x in reference))
    exact = [(i,j) for i,j in pairs if observed[i][2] == reference[j][1]]
    identity = len(exact)/len(pairs) if pairs else 0
    coverage = len(exact)/len(observed)
    report = {
        "source": "GPCRdb REST residues/extended", "entry": entry, "chain": chain,
        "residue_count": len(observed), "aligned": len(pairs), "exact_matches": len(exact),
        "identity": round(identity, 4), "coverage": round(coverage, 4),
        "mismatches": [
            {"pdb_residue": observed[i][0]+observed[i][1],
             "pdb_aa": observed[i][2], "reference_position": reference[j][0],
             "reference_aa": reference[j][1]} for i,j in pairs if observed[i][2] != reference[j][1]
        ]
    }
    pathlib.Path(audit_path).write_text(json.dumps(report, indent=2))
    if len(exact) < 30 or identity < .85 or coverage < .65:
        raise ValueError(f"Alignment rejected (identity={identity:.1%}, coverage={coverage:.1%}); audit: {audit_path}")
    if any(observed[i][1] for i,j in exact):
        raise ValueError("PDB insertion codes present. Numeric-only LigInt mapping cannot represent them safely.")
    rows = ["chain\tresi\taa\tgpcrdb\tsegment"]
    for i,j in exact:
        position, aa, generic, segment = reference[j]
        if generic:
            rows.append(f"{chain}\t{observed[i][0]}\t{aa}\t{generic}\t{segment}")
    if len(rows) <= 10:
        raise ValueError("Too few GPCRdb-numbered positions to write")
    pathlib.Path(output).write_text("\n".join(rows) + "\n")
    report["numbered_positions"] = len(rows)-1
    pathlib.Path(audit_path).write_text(json.dumps(report, indent=2))
    return report

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", default="auto", help="PDB chain ID, or auto to identify the best matching receptor chain")
    ap.add_argument("--entry", default="auto", help="GPCRdb entry or auto (infer from PDB identifier)")
    ap.add_argument("--pdb-id", help="Four-character PDB ID when --entry auto")
    ap.add_argument("--output", required=True)
    ap.add_argument("--audit")
    ap.add_argument("--cache-dir", default=str(pathlib.Path.home()/".cache"/"ligint"/"gpcrdb"))
    a = ap.parse_args()
    if not all(ch.isalnum() or ch == "_" for ch in a.entry):
        ap.error("Invalid GPCRdb entry")
    audit = a.audit or a.output + ".audit.json"
    try:
        if a.entry == "auto":
            a.entry = identify_entry_from_pdb(a.pdb_id or pathlib.Path(a.pdb).stem[:4])
        report = create(a.pdb, a.chain, a.entry, a.cache_dir, a.output, audit)
    except Exception as exc:
        print(f"GPCRdb mapping failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2))
    return 0

if __name__ == "__main__":
    sys.exit(main())
