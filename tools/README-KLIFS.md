# KLIFS -> LigInt pocket definition (experimental integration)

LigInt v0.7.0 already accepts `--pocket FILE` and renders CSV/JSON/SVG fingerprints.
This helper obtains an **85-position kinase pocket mapping** from KLIFS and exports
its available PDB residues using exactly that existing file format.

## Online example

Install `opencadd` with its compatible dependencies into a Python environment.
Then run:

```bash
python tools/klifs_to_pocket.py --pdb 3EYG --chain A --model B --output 3EYG-klifs-pocket.tsv
ligint 3EYG.pdb --fingerprint svg --pocket 3EYG-klifs-pocket.tsv > 3EYG.svg
```

The first step needs internet access. `--model` refers to the **KLIFS alternate
model**, which is not necessarily a PDB MODEL record. To see matching choices,
omit `--model` if ambiguous; the tool will report available mappings.

## Offline import of KLIFS mapping

```bash
python tools/klifs_to_pocket.py --pdb 3EYG --chain A \
  --mapping-csv klifs-pocket-export.csv --output 3EYG-klifs-pocket.tsv
```

CSV/TSV columns: `residue.klifs_id`, `residue.id`, `residue.klifs_region`
(optionally `residue.aa`). IDs must refer to the **PDB author numbering**
of the same coordinate file and chain. Do not substitute mmCIF `label_seq_id`
for PDB author identifiers without mapping. Missing positions are omitted,
not filled with fabricated residue numbers. Validate chain and numbering in
LigInt's structure before relying on this output.

This produces **structural annotations**, not KLIFS-computed ligand contacts.
LigInt calculates contacts independently from the structure file. Use the same
crystallographic model, chain and ligand instance for both. Where KLIFS contains
multiple chains or alternate models, select the matching one explicitly.

Caveat: The online implementation uses the documented OpenCADD KLIFS client API,
but online service integration was not end-to-end tested in this isolated build
sandbox. Offline conversion and TSV compatibility can be tested locally.
