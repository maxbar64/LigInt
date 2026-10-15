# Quick start

This page takes you from a structure file to a table of interactions in a few minutes. The example is the structure
**1AZM** (human carbonic anhydrase II with the inhibitor acetazolamide, ligand name `AZM`, and a zinc ion), which you can
download from the [RCSB Protein Data Bank](https://www.rcsb.org/structure/1AZM) as `1AZM.pdb` or `1AZM.cif`.

## In the graphical interface

1. Start **ligint** from the applications menu (Linux), the Start menu (Windows) or with `ligint-gui` in a terminal.
2. Click **Open structure...** (or press `Ctrl+O`, or drag the file onto the window) and choose `1AZM.pdb`.
3. The **Ligand** list now shows `AZM A 262 (13 atoms)`. Ions and very small molecules (fewer than 6 atoms) are hidden unless you
   tick **Also show ions and small molecules**.
4. Press **Analyze**. The ligand and the residues that touch it appear in the 3D view; the **Interactions** tab lists each
   interaction with its atoms and distance.
5. Drag in the view to rotate, use the wheel to zoom, click a row of the table to highlight that interaction.
6. **File > Export interactions (CSV)...** saves the table, **File > Save image...** saves the view as a PNG, and
   **File > Copy command line** copies the `ligint` command that reproduces the analysis.

![The graphical interface on 1HSG](img/gui-1HSG-window.png)

The full description of the window is in [The graphical interface](gui.md).

## On the command line

Without a ligand name, `ligint` lists the ligands of the file:

```text
$ ligint 1AZM.pdb
Ligands (HETATM) present: AZM(1) ZN(1)
```

Give the name of a ligand to analyze it:

```text
$ ligint 1AZM.pdb AZM
Ligand annotation AZM 262A: geometric heuristics
Mean B-factor: ligand 8.6, rest of the structure 15.4

=== Ligand AZM 262 chain A: 13 atoms, 1 rings ===
Type                     Fav Ligand atom          Protein atom         d(A)   Notes
Hydrogen bond            +   AZM 262A:N1          THR 199A:OG1         3.00   
Hydrogen bond            +   AZM 262A:O2          THR 199A:N           3.08   
Hydrophobic contact      +   AZM 262A:C4          PHE 91A:CZ           3.78   
Hydrophobic contact      +   AZM 262A:C4          PHE 91A:CE2          3.82   
Metal coordination       +   AZM 262A:N1          ZN 261A:ZN           2.01   
Water bridge             +   AZM 262A:N4          GLN 92A:OE1          2.99   via HOH 360A:O; d(H2O-prot) 3.32

Summary by type:
  Hydrogen bond                2
  Hydrophobic contact          2
  Metal coordination           1
  Water bridge                 1

Summary by protein residue:
  PHE 91A       2  Hydrophobic contact
  THR 199A      2  Hydrogen bond
  GLN 92A       1  Water bridge
  ZN 261A       1  Metal coordination

Coordination sphere of the metals in contact with the ligand:
  ZN 261A:ZN 4 atoms: AZM 262A:N1 (2.01) HIS 119A:ND1 (2.02) HIS 94A:NE2 (2.04) HIS 96A:NE2 (2.17)
```

How to read it:

- Each row of the table is one interaction: its **type**, whether it is favorable (`+`) or not (`-`, clashes and repulsions),
  the **ligand atom**, the **protein atom** and the **distance** in angstroms, and **notes** (angles, the water that mediates a bridge...).
- The **summaries** count the interactions by type and by protein residue.
- The last block lists the **coordination sphere** of the metals that the ligand touches. Here the ligand nitrogen and three
  histidines hold the zinc.

Useful variations:

```text
ligint 1AZM.pdb AZM --csv > interactions.csv          the same table as CSV, for a spreadsheet or a script
ligint 1AZM.pdb AZM A 262                             only chain A, residue 262 (when the ligand appears several times)
ligint 1AZM.pdb AZM --no-water                        ignore water molecules (no water bridges)
ligint 1AZM.pdb AZM --fingerprint svg > figure.svg    the interaction fingerprint as a figure
ligint 1AZM.cif AZM                                   the same analysis from the mmCIF file
```

All the options are in the [Command-line reference](cli.md).

## Where to get structures

Search by name or code at [rcsb.org](https://www.rcsb.org) and download the **PDB** or **PDBx/mmCIF** file. Files from the archive
are often compressed (`.gz`): decompress them first (`gunzip 1abc.cif.gz`). New entries of the Protein Data Bank will only be
distributed as PDBx/mmCIF: see [Input formats](formats.md).

## Things to know before trusting a result

- Crystal structures do not contain hydrogens, so the donor and acceptor roles and the charges of the ligand are **inferred from the
  geometry** (for example, carboxylates are treated as anionic and aliphatic amines as cationic). An SDF or MOL2 file of the ligand
  can supply the exact chemistry: see [Interactions and assumptions](interactions.md).
- The results are **geometric criteria applied to a model**, not a measurement of binding. [Validation and limits](validation.md)
  explains how they compare with another tool and where they are weak.
