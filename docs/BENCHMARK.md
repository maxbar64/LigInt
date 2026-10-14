# Comparison with PLIP

This document describes how the agreement between `ligint` and [PLIP](https://github.com/pharmai/plip) was measured
(version 3.0.1; on the 2.4.0 report of 1AZM provided by the user it gives the same interactions), what can and cannot be
concluded, and how to repeat the measurement.

**PLIP is not the truth.** The comparison measures the agreement with another tool, which has its own conventions
(protonation with Open Babel, different thresholds, filters on water bridges). A disagreement is a reason to look at
the case, not an automatic error. There is no comparison with independent experimental data.

## Method

- **Structures:** 150 PDB complexes in four groups (A and B with 25 each, C and D with 50 each; list in `groups.tsv`),
  chosen to cover different protein families and chemistries (kinases, proteases, nuclear receptors, GPCRs, enzymes with
  cofactors and heme, metalloproteins, macrocycles, sugars, halogenated ligands...). The IDs were chosen from memory and
  checked against the title after the download.
- **Sites:** for each PDB, PLIP analyzes all the small-molecule ligands; ligint is run on the same ligand (chain,
  residue). Only **one copy per ligand and PDB** is counted (otherwise structures with 20 copies dominate the totals).
  Ligands classified by PLIP as polymers (long glycan chains), peptides and nucleic acids are excluded.
- **Comparison at the level of residue and type:** (type, residue, number, chain). Types compared: hydrophobic, H-bond,
  water bridge, salt bridge, π-stacking, cation-π, halogen. Not compared: clashes, repulsions, covalent bonds, contacts
  with the backbone, metal coordination (PLIP describes it between protein and metal, ligint between ligand and metal).
- **Convention for H-bonds between opposite charges:** PLIP reports a salt bridge and an H-bond on the same pair together
  in only 13% of cases; the "charge-assisted" H-bonds of ligint (or those merged into the salt bridge) do not enter the
  H-bond score.
- **Sensitivity** = PLIP interactions found by ligint; **precision** = ligint interactions also present in PLIP. The 95%
  intervals are bootstrap over the structures (the interactions of one structure are not independent).

## Groups and use during development

| Group | Role |
|---|---|
| A, B | used to look for errors and correct the program (heme, water bridges, linked residues, amines, amidines...) |
| C | first group not used for the corrections; its disagreements were nevertheless examined |
| D | not used for any decision except the check of the options `--hb-long` and `--salt-centroid` (proposed after seeing A, B and C) |

The results for A, B and C are therefore optimistic compared with unseen complexes; the most honest estimate is D.

## Results (default thresholds, one copy per ligand)

| Group | Sensitivity | Precision |
|---|---|---|
| A | 80% [73–87] | 72% [68–77] |
| B | 82% [78–85] | 76% [72–81] |
| C | 77% [73–80] | 70% [66–73] |
| D | 75% [70–80] | 70% [66–74] |
| D with `--hb-long 4.1 --salt-centroid 5.5` | 79% [74–84] | 70% [65–74] |

Group D by type (default):

| Type | Sensitivity | Precision |
|---|---:|---:|
| hydrophobic | 86% | 81% |
| hydrogen bond | 72% | 75% |
| water bridge | 65% | 40% |
| salt bridge | 52% | 81% |
| π-stacking | 94% | 83% |
| cation-π | 67% | 44% |
| halogen | 40% | 50% |

Effect of the two options on D (paired difference, bootstrap over the structures, options − default):

| | Sensitivity | Precision |
|---|---|---|
| H-bond (`--hb-long 4.1`) | +9.9 points [+5.9; +14.8] | −2.7 [−5.0; −0.3] |
| Salt bridges (`--salt-centroid 5.5`) | +5.1 [0.0; +11.8] | −2.3 [−9.3; +2.5] |
| Total | +4.2 [+2.7; +5.9] | −0.8 [−1.8; +0.1] |

The options are therefore off by default: `--hb-long` has a clear and confirmed effect, `--salt-centroid` an uncertain one.

## Why ligint and PLIP diverge

- **Fluorine:** PLIP counts fluorine as a halogen-bond donor too (12 of the 32 cases of the first two benchmarks); ligint only Cl, Br, I.
- **Histidine:** PLIP considers it positively charged; about 36% of the salt bridges missed by ligint involve it.
- **Long H-bonds:** about two thirds of the contacts found only by PLIP are between 3.5 and 4.1 Å (the PLIP limit, with
  computed H geometry and a much looser angular constraint); `--hb-long` recovers some of them.
- **Hydrophobic:** PLIP counts carbons of polar residues (Gln, Asp, Glu, Arg) and of heme; ligint limits itself to hydrophobic residues.
- **Water bridges:** PLIP applies its own filters (a donor already in a direct H-bond does not form a bridge; at most two
  bridges per donor water); ligint does not, and reports many more bridges (precision ≈ 40%).
- **Ligand protonation:** ligint assumes aliphatic amines protonated and carboxylates anionic; PLIP uses Open Babel.

## Repeating the measurement

```
pip install plip openbabel-wheel
python3 tests/benchmark.py folder_with_pdb/ --groups groups.txt --out results/ [--all-copies] ["--extra=--plip-like"]
```

It writes `report.md`, `disagreements.csv` (every disagreement with atoms, distances and an automatic probable cause)
and `per_category.csv` (TP/FP/FN per site and type, for bootstrap and paired comparisons). The reports in this
repository are in `docs/benchmark/`: `AB`, `C`, `D_default`, `D_hb-long_salt-centroid`. The probable causes are an
automatic classification of convenience (distance beyond the threshold, histidine partner...).

## Limits

- A single reference tool; no experimental truth.
- No comparison for DNA/RNA (PLIP ignores them), halogens with few cases (order of 10), cation-π few cases.
- The IDs were chosen from memory by a single person: the representativeness of the groups was not verified.
