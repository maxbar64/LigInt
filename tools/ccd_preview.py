#!/usr/bin/env python3
"""Fetch CCD bond topology for ligand IDs, cache CIF, emit atom-name keyed JSON.
Uses only stdlib; read-only remote source is RCSB CCD individual CIF.
"""
import argparse
import json
from pathlib import Path
import re
import shlex
from urllib.request import Request, urlopen

ORDER={"SING":1,"DOUB":2,"TRIP":3,"AROM":4}

def loop_records(lines, start, columns):
    """Read CIF loop values across physical lines, rather than assuming one row per line."""
    values = []
    i = start
    while i < len(lines):
        raw = lines[i].strip()
        if raw == "#":
            i += 1
            break
        if raw in ("loop_", "stop_") or raw.startswith(("data_", "save_", "_")):
            break
        if raw.startswith(";"):
            raise ValueError("Multiline CIF text fields in CCD atom/bond loop are unsupported")
        if raw:
            values.extend(shlex.split(raw, comments=False, posix=True))
        i += 1
    if len(values) % len(columns):
        raise ValueError("Incomplete CCD loop row (field count mismatch)")
    return (values[j:j + len(columns)] for j in range(0, len(values), len(columns)))


def parse_component(text, component):
    lines=text.splitlines()
    result=[]
    i=0
    while i<len(lines):
        if lines[i].strip()!="loop_":
            i+=1;continue
        i+=1
        columns=[]
        while i<len(lines) and lines[i].lstrip().startswith("_"):
            columns.append(lines[i].strip().split()[0])
            i+=1
        if not columns or not any(c.startswith("_chem_comp_bond.") for c in columns):
            continue
        target=["_chem_comp_bond.atom_id_1","_chem_comp_bond.atom_id_2",
                "_chem_comp_bond.value_order"]
        if not all(x in columns for x in target):
            continue
        aromatic="_chem_comp_bond.pdbx_aromatic_flag"
        indices=[columns.index(x) for x in target]
        arom_index=columns.index(aromatic) if aromatic in columns else None
        comp_col=columns.index("_chem_comp_bond.comp_id") if "_chem_comp_bond.comp_id" in columns else None
        for values in loop_records(lines, i, columns):
            if comp_col is None or values[comp_col].upper()==component:
                order=ORDER.get(values[indices[2]].upper())
                if order is None:raise ValueError("Unsupported bond order")
                if arom_index is not None and values[arom_index].upper()=="Y":order=4
                result.append([values[indices[0]],values[indices[1]],order])
        break
    if not result:raise ValueError("No CCD bond table found")
    return result

def parse_atom_elements(text):
    """Return atom-name -> element from the CCD atom category; fail closed."""
    lines=text.splitlines()
    i=0
    while i<len(lines):
        if lines[i].strip()!="loop_":
            i+=1;continue
        i+=1;headers=[]
        while i<len(lines) and lines[i].lstrip().startswith("_"):
            headers.append(lines[i].strip().split()[0]);i+=1
        names=("_chem_comp_atom.atom_id","_chem_comp_atom.type_symbol")
        if not all(k in headers for k in names):continue
        name_idx= headers.index(names[0]);elem_idx=headers.index(names[1])
        result={}
        for fields in loop_records(lines, i, headers):
            if fields[name_idx] in result:raise ValueError("Duplicate CCD atom")
            result[fields[name_idx]]=fields[elem_idx].upper()
        if result:return result
    raise ValueError("CCD atom definitions missing")

def donor_hydrogens(bonds, elements):
    """CCD-defined N-H/O-H sites; CCD is an explicit *component* model, not experimental pH."""
    sites={}
    for a,b,_ in bonds:
        if elements.get(a)=="H" and elements.get(b) in ("N","O"):
            sites[b]=sites.get(b,0)+1
        elif elements.get(b)=="H" and elements.get(a) in ("N","O"):
            sites[a]=sites.get(a,0)+1
    return sites

def get_component(component, cache):
    if not re.fullmatch(r"[A-Z0-9]{1,8}",component):raise ValueError("Invalid CCD ID")
    path=cache/(component+".cif")
    if path.exists():
        text=path.read_text()
    else:
        url=f"https://files.rcsb.org/ligands/download/{component}.cif"
        with urlopen(Request(url,headers={"User-Agent":"LigInt-CCD-preview/1.0"}),timeout=15) as response:
            raw=response.read(2_000_001)
        if len(raw)>2_000_000:raise ValueError("CCD entry exceeds size limit")
        text=raw.decode("utf-8")
    bonds=parse_component(text,component)
    if not path.exists():path.write_text(text)
    elements=parse_atom_elements(text)
    return bonds,donor_hydrogens(bonds,elements),elements

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("components",nargs="+")
    p.add_argument("--cache-dir",default=str(Path.home()/".cache"/"ligint"/"ccd"))
    args=p.parse_args()
    cache=Path(args.cache_dir)
    cache.mkdir(parents=True,exist_ok=True)
    output={}
    for component in dict.fromkeys(s.upper() for s in args.components):
        try:
            bonds,hydrogens,elements=get_component(component,cache)
            output[component]={"bonds":bonds,"donor_hydrogens":hydrogens,"elements":elements}
        except Exception as exc:output[component]={"error":str(exc)}
    print(json.dumps(output))
if __name__=="__main__":main()
