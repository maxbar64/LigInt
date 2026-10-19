#!/usr/bin/env python3
"""Non-destructive ligand-blind hierarchical grouping of predicted pocket voxels."""
import argparse
import json
import math
from collections import defaultdict
from pocket_contact_audit import read_pockets, grid_adjacency


def candidates(pockets, spacing, min_faces=8, min_fraction=0.08,
               max_volume=1800, max_members=4):
    faces = grid_adjacency(pockets, spacing)
    sizes = {p: len(v) for p, v in pockets.items()}
    parent = {p:p for p in pockets}
    groups = {p:{p} for p in pockets}
    history = []

    def root(p):
        while parent[p]!=p:
            p=parent[p]
        return p

    edges = []
    for (a,b),n in faces.items():
        # Contact ratio normalized to the smaller subpocket's six-neighbor face budget.
        ratio=n/(6*min(sizes[a],sizes[b]))
        if n>=min_faces and ratio>=min_fraction:
            edges.append((-ratio,-n,a,b))
    edges.sort()
    for neg_ratio,neg_n,a,b in edges:
        ra,rb=root(a),root(b)
        if ra==rb:
            continue
        left,right=groups[ra],groups[rb]
        merged_volume=sum(sizes[p] for p in left|right)*spacing**3
        if merged_volume>max_volume or len(left|right)>max_members:
            continue
        representative=min(ra,rb)
        other=max(ra,rb)
        parent[other]=representative
        groups[representative]=left|right
        del groups[other]
        history.append({"joined_edge":[a,b],"contact_faces":-neg_n,
                        "normalized_contact":round(-neg_ratio,5),
                        "members":sorted(groups[representative]),
                        "volume_a3":round(merged_volume,3)})
    return {"settings":{"spacing_a":spacing,"min_contact_faces":min_faces,
                        "min_contact_fraction":min_fraction,"max_group_volume_a3":max_volume,
                        "max_group_members":max_members},
            "subpocket_count":len(pockets),"groups":[{"members":sorted(m),
                 "volume_a3":round(sum(sizes[p] for p in m)*spacing**3,3)}
                 for m in sorted(groups.values(),key=lambda m:(-sum(sizes[p] for p in m),min(m)))
                 if len(m)>1],"merge_history":history,
            "note":"Ligand-blind candidate grouping; original pockets unchanged. Validate with independent crystals."}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pockets-pdb",required=True)
    ap.add_argument("--spacing",type=float,default=1)
    ap.add_argument("--min-contact-faces",type=int,default=8)
    ap.add_argument("--min-contact-fraction",type=float,default=.08)
    ap.add_argument("--max-group-volume",type=float,default=1800)
    ap.add_argument("--max-group-members",type=int,default=4)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()
    if not math.isfinite(a.spacing) or a.spacing<=0 or not math.isfinite(a.min_contact_fraction) or not 0<=a.min_contact_fraction<=1 or not math.isfinite(a.max_group_volume) or a.max_group_volume<=0 or a.min_contact_faces<1 or a.max_group_members<2:
        ap.error("invalid grouping parameters")
    result=candidates(read_pockets(a.pockets_pdb),a.spacing,a.min_contact_faces,a.min_contact_fraction,a.max_group_volume,a.max_group_members)
    with open(a.output,"w",encoding="utf-8") as f:
        json.dump(result,f,indent=2)
        f.write("\n")


if __name__=="__main__":
    main()
