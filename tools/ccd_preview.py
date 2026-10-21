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
        while i<len(lines):
            line=lines[i].strip()
            if not line or line=="#":
                i+=1
                if line=="#":break
                continue
            if line=="loop_" or line.startswith("_") or line.startswith("data_") or line.startswith("save_"):
                break
            if line.startswith(";"):
                raise ValueError("Unsupported multiline bond entry")
            values=shlex.split(line,comments=False,posix=True)
            if len(values)!=len(columns):
                raise ValueError("CCD bond row has unexpected number of fields")
            if comp_col is None or values[comp_col]==component:
                order=ORDER.get(values[indices[2]].upper())
                if order is None:raise ValueError("Unsupported bond order")
                if arom_index is not None and values[arom_index].upper()=="Y":order=4
                result.append([values[indices[0]],values[indices[1]],order])
            i+=1
    if not result:raise ValueError("No CCD bond table found")
    return result

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
    return bonds

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("components",nargs="+")
    p.add_argument("--cache-dir",default=str(Path.home()/".cache"/"ligint"/"ccd"))
    args=p.parse_args()
    cache=Path(args.cache_dir)
    cache.mkdir(parents=True,exist_ok=True)
    output={}
    for component in dict.fromkeys(s.upper() for s in args.components):
        try:output[component]={"bonds":get_component(component,cache)}
        except Exception as exc:output[component]={"error":str(exc)}
    print(json.dumps(output))
if __name__=="__main__":main()
