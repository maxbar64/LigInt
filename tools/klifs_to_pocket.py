#!/usr/bin/env python3
"""Create LigInt --pocket TSV from a KLIFS 85-position mapping.

Online mode requires opencadd (with its KLIFS client); offline mode reads a CSV
with columns residue.klifs_id, residue.id, residue.klifs_region, etc.
The output is a LigInt pocket *definition*, not a calculated interaction set.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

COLORS = {
    'I':'#e7dc71', 'g.l':'#8acb7b', 'II':'#f4cf65', 'III':'#e6ab63',
    'aC':'#e47875', 'b.l':'#7bc96f', 'IV':'#eee679', 'V':'#f0e877',
    'GK':'#efcb64', 'hinge':'#cc83d8', 'linker':'#80dbe1',
    'aD':'#df7f82', 'aE':'#d77a75', 'aF':'#da7778',
    'VI':'#f4e97d', 'VII':'#f3e579', 'VIII':'#f3e579',
    'c.l':'#edb164', 'DFG':'#8ea5ec', 'xDFG':'#8ea5ec', 'a.l':'#8ea5ec',
}
AA3 = dict(zip('ACDEFGHIKLMNPQRSTVWY',
    'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split()))

def clean_region(region):
    # Keep KLIFS-specific labels, which are also a valid LigInt TSV field.
    return str(region).strip().replace(' ', '_')

def csv_rows(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        first = f.readline()
        f.seek(0)
        return list(csv.DictReader(f, delimiter='\t' if '\t' in first else ','))

def pick(record, *fields):
    for f in fields:
        v = record.get(f)
        if v is not None and str(v).strip() not in ('', '_', 'None', 'nan', 'NaN'):
            return str(v).strip()
    return ''

def build(rows, chain, name, pdb_id, strict=True):
    out = [f'# KLIFS-derived pocket definition for {pdb_id.upper()} chain {chain}',
           '# Columns: chain  PDB author residue id  KLIFS position  KLIFS region  optional residue name',
           '# Positions missing in the experimental structure are omitted; no guessed numbering.',
           f'@name {name} (KLIFS, {pdb_id.upper()} chain {chain})']
    parsed = []
    seen = set()
    for row in rows:
        kid = pick(row, 'residue.klifs_id', 'klifs_id', 'klifs_position', 'position')
        pdbres = pick(row, 'residue.id', 'residue_id', 'pdb_residue', 'resi')
        region = pick(row, 'residue.klifs_region', 'klifs_region', 'region')
        aa = pick(row, 'residue.aa', 'residue.amino_acid', 'residue_name', 'resname')
        if not kid.isdigit() or int(kid) not in range(1, 86):
            if strict: raise ValueError(f'Invalid KLIFS position: {kid!r}')
            continue
        if kid in seen: raise ValueError(f'Duplicate KLIFS position {kid}')
        seen.add(kid)
        if not pdbres: continue
        if not re.fullmatch(r'-?\d+[A-Za-z]?', pdbres):
            raise ValueError(f'Unusable PDB residue identifier {pdbres!r}; must be author residue number')
        if not region: raise ValueError(f'Missing KLIFS region for position {kid}')
        if len(aa) == 1: aa = AA3.get(aa.upper(), '')
        if aa and not re.fullmatch(r'[A-Za-z0-9]{3}', aa): aa = ''
        parsed.append((int(kid), pdbres, clean_region(region), aa.upper()))
    if not parsed: raise ValueError('No mapped KLIFS residues found')
    if len(parsed) > 85: raise ValueError('Unexpected >85 mapped positions')
    for region in dict.fromkeys(r[2] for r in sorted(parsed)):
        color = COLORS.get(region, '#bdc6d0')
        out.append(f'@color {region} {color}')
    for kid, pdbres, region, aa in sorted(parsed):
        out.append('\t'.join([chain, pdbres, str(kid), region] + ([aa] if aa else [])))
    return '\n'.join(out) + '\n', len(parsed)

def online(pdb, chain, model):
    try:
        from opencadd.databases.klifs import setup_remote
    except ImportError as exc:
        raise RuntimeError('Online mode requires opencadd; install a compatible release or use --mapping-csv') from exc
    remote = setup_remote()
    structures = remote.structures.by_structure_pdb_id([pdb.lower()])
    if structures is None or structures.empty:
        raise ValueError(f'{pdb} not present in KLIFS')
    # KLIFS can distinguish chains and alternate models within a single PDB code.
    rows = structures.loc[structures['structure.chain'].astype(str).str.upper() == chain.upper()]
    if model is not None:
        rows = rows.loc[rows['structure.alternate_model'].astype(str).str.upper() == model.upper()]
    if len(rows) != 1:
        choices = structures[['structure.klifs_id', 'structure.chain', 'structure.alternate_model']].to_dict('records')
        raise ValueError(f'Expected one KLIFS model for {pdb} chain {chain} model {model}; matches={len(rows)}. Available: {choices}')
    structure_id = int(rows.iloc[0]['structure.klifs_id'])
    pocket = remote.pockets.by_structure_klifs_id(structure_id)
    if pocket is None or pocket.empty:
        raise ValueError(f'KLIFS pocket mapping empty for structure {structure_id}')
    # Do not fabricate residue identities from the consensus KLIFS pocket sequence.
    return pocket.to_dict('records'), structure_id

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pdb', required=True, help='four-character PDB ID, e.g. 3EYG')
    p.add_argument('--chain', required=True, help='author chain ID, e.g. A')
    p.add_argument('--model', help='KLIFS alternate model, e.g. B; required if multiple KLIFS matches')
    p.add_argument('--mapping-csv', type=Path, help='offline: KLIFS pocket export CSV/TSV')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    try:
        if not re.fullmatch(r'[A-Za-z0-9]{4}', args.pdb): raise ValueError('PDB ID must be four alphanumeric characters')
        if not re.fullmatch(r'[A-Za-z0-9]+', args.chain): raise ValueError('Invalid chain ID')
        if args.mapping_csv:
            rows, source = csv_rows(args.mapping_csv), str(args.mapping_csv)
        else:
            rows, source = online(args.pdb, args.chain, args.model)
        tsv, n = build(rows, args.chain, 'Kinase pocket', args.pdb)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(tsv, encoding='utf-8')
        print(f'Wrote {n}/85 mapped KLIFS positions to {args.output} (source: {source})')
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        p.exit(2, f'KLIFS pocket generation failed: {exc}\n')

if __name__ == '__main__': main()
