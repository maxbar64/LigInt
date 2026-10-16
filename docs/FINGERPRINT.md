# Interaction fingerprints

A fingerprint summarizes the interactions of one ligand as a table with **one row per protein residue and one column
per interaction class**. It can be exported to CSV and JSON, drawn as an SVG figure (residues grouped by pocket region,
one glyph per class) and shown in the GUI. It is the kind of table/figure used to compare how different ligands use the
same binding pocket.

![Fingerprint of imatinib in 1IEP, without a pocket definition](img/fingerprint-1IEP.svg)

## Classes

| Column | Meaning | Derived from |
|---|---|---|
| `hydrophobic` | hydrophobic contact | *Hydrophobic contact* (and *Apolar contact (backbone)*) |
| `face_to_face` | aromatic face-to-face | *pi-stacking (parallel)* |
| `edge_to_face` | aromatic edge-to-face | *pi-stacking (T-shaped)* |
| `hbond_donor` | the **protein residue donates** the hydrogen | *Hydrogen bond*, protein atom is the donor |
| `hbond_acceptor` | the **protein residue accepts** the hydrogen | *Hydrogen bond*, protein atom is the acceptor |
| `ionic_positive` | salt bridge with a positively charged residue | *Salt bridge*, charge of the protein atom |
| `ionic_negative` | salt bridge with a negatively charged residue | idem |
| `cation_pi` | cation-pi | *Cation-pi* |
| `halogen` | halogen bond | *Halogen bond* |
| `metal` | metal coordination | *Metal coordination* |
| `water_bridge` | water-mediated contact | *Water bridge* |

The first seven classes are the ones of the usual kinase-ligand fingerprint figures; the others are extras. Each cell of
the table is a **count of contacts** (the JSON also gives the shortest distance per class). A salt bridge that is also an
H-bond counts in the ionic class and in the donor/acceptor class. Clashes, repulsions and possible covalent bonds are
not part of the fingerprint.

**Caveats.** The donor/acceptor role and the charge of each residue come from the same geometric heuristics as the
interactions (histidine is never charged, amines and carboxylates follow the assumptions in the README). The classes are
*not* guaranteed to coincide with the definitions of other tools: measure before treating them as equivalent.

## Command line

```
ligint structure.pdb LIG --fingerprint csv                       one row per residue in contact
ligint structure.pdb LIG --fingerprint json
ligint structure.pdb LIG --fingerprint svg > figure.svg
ligint structure.pdb LIG --fingerprint svg --pocket pocket.tsv   with a pocket definition
```

- `--fingerprint` replaces the interaction table (the service information goes to stderr, as with `--csv`).
- CSV and JSON contain every ligand found (a `ligand` column / one object per ligand). The SVG draws only the first
  ligand and warns on stderr if there are more: select one with the chain and residue number.
- `--pocket` needs `--fingerprint`. Unknown formats, missing files and `--pocket` alone end with an error (code 1);
  malformed lines of the pocket file are skipped with a warning.

CSV columns: `ligand,residue,chain,resi,position,region` followed by one count column per class
(`hydrophobic,face_to_face,edge_to_face,hbond_donor,hbond_acceptor,ionic_positive,ionic_negative,cation_pi,halogen,metal,water_bridge`).

## Pocket definition file

The pocket definition is **data, not code**: a text file that says which residues form the pocket, in which order and
how they are grouped. Without it the fingerprint lists only the residues in contact, sorted by chain and number. With it,
every position is listed (also those without interactions, as in the usual figures) in the order of the file, and the
contact residues that are not in the definition are appended and flagged (`in_pocket` = false).

```
# comment
@name HIV-1 protease pocket (example)
@color catalytic #e5857f
A	25	A25	catalytic	ASP
A	50	A50	flap	ILE
```

- One position per line: **chain, residue number, position label, region [, residue name]**. Fields are separated by
  tabs or commas (then they may contain spaces) or, if neither is present, by whitespace.
- The residue number accepts an insertion code (`134A`). The chain `*` (or `-`, `_`) matches any chain.
- The region groups consecutive positions in the figure (a colored bar); the order of the file is the drawing order.
- The residue name is optional: it labels the positions that have no contact, and a mismatch with the structure is
  reported as a warning.
- `@name` sets the title and `@color REGION #rrggbb` sets the color of a region; other `@` directives are ignored with a warning.
- A numbering is specific to a structure: a different structure of the same protein needs its own file (or an
  alignment-based mapping, which is not part of the tool yet).

An illustrative file is in [`examples/hiv-protease-pocket.tsv`](examples/hiv-protease-pocket.tsv) (the regions and
labels are an example of the format, **not** a reference mapping). It produces:

![Fingerprint of indinavir in 1HSG with the example definition](img/fingerprint-1HSG.svg)

## GUI

The **Fingerprint** tab shows the figure and the table. *Load definition...* in the control panel loads a pocket file
(*Clear* removes it); the tab updates immediately. **File > Export fingerprint...** saves SVG, CSV or JSON (chosen by the
extension). The figure needs the Qt Svg module; without it the GUI shows the table only.

![GUI, Fingerprint tab](img/gui-1HSG-fingerprint.png)

## Library

```cpp
#include <ligint/fingerprint.hpp>

ligint::Result r = ligint::analyzeFile("1HSG.pdb", opt);
ligint::PocketDefinition pocket;
std::vector<std::string> warnings;
ligint::loadPocketDefinition("pocket.tsv", pocket, warnings);          // optional
ligint::Fingerprint fp = ligint::makeFingerprint(r.ligands[0], &pocket);
for (const auto& e : fp.entries)
    std::cout << e.residue << " " << e.position << " " << e.counts[ligint::FpHBondDonor] << "\n";
std::string svg = ligint::fingerprintSvg(fp);
```

`writeFingerprintCsv` and `writeFingerprintJson` take a vector of fingerprints (one per ligand). The library does not
depend on Qt.

## Not covered yet

- Mapping residues to pocket positions: the position of each residue in a reference pocket (as defined for kinases, for
  example) is not computed; it must come from the definition file.
- Comparison of the fingerprint with those of other tools.

The fingerprints of several ligands on the same columns, one bit for each position and class, are in the [interaction matrix](MATRIX.md).
