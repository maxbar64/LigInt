# ligint manual

**ligint** finds the non-covalent interactions between a small-molecule **ligand** and the **protein** around it in a structure file
(PDB or PDBx/mmCIF): hydrogen bonds, salt bridges, hydrophobic contacts, pi-stacking, cation-pi, halogen bonds, metal coordination and
water bridges, plus steric clashes and electrostatic repulsions. It can also summarize them as an **interaction fingerprint**: a table
with one row per protein residue and one column per interaction class, drawn as a figure.

It comes as a **command-line program** (`ligint`) and a **graphical interface** (`ligint-gui`) that share the same analysis engine.
The graphical interface shows the ligand and its surroundings in 3D, lists the interactions in tables, and writes, for every
analysis, the command line that reproduces it outside the program.

## Contents

- [Quick start](quickstart.md): your first analysis, in the graphical interface and on the command line.
- [The graphical interface](gui.md): every control, view and export.
- [Command-line reference](cli.md): arguments, options, thresholds, output formats, exit codes, examples.
- [Interactions and assumptions](interactions.md): what is recognized, with which criteria, and the chemical assumptions behind them.
- [Input formats](formats.md): PDB and PDBx/mmCIF, alternate conformations, chains of several characters.
- [Interaction fingerprints](fingerprint.md): classes, pocket definition files, CSV, JSON and SVG output.
- [Interaction matrix](matrix.md): the fingerprints of several ligands on the same columns, as CSV, JSON or a heat map.
- [KLIFS kinase pockets](klifs.md): pocket definitions for kinase fingerprints, from a KLIFS mapping.
- [GPCR pockets (GPCRdb)](gpcrdb.md): pocket definitions for GPCRs, from the generic residue numbers of GPCRdb.
- [Downloading from the RCSB](rcsb.md): fetch a structure by its PDB ID, in the interface and on the command line.
- [Validation and limits](validation.md): how the results were checked, and what to keep in mind when reading them.
- [Comparison with PLIP](benchmark.md): the measurement of the agreement with another tool.
- [Installation](install.md): packages and installers for Linux and Windows.
- [Troubleshooting](troubleshooting.md): common problems and their causes.
- [For developers](developers.md): architecture, tests, how to add an interaction type.
- [Changelog](changelog.md) and [About](about.md).

## Where to start

| If you want to... | Go to |
|---|---|
| analyze a structure by clicking | [Quick start](quickstart.md), then [The graphical interface](gui.md) |
| script many analyses | [Command-line reference](cli.md) |
| know exactly what counts as an interaction | [Interactions and assumptions](interactions.md) |
| make a pocket-position figure for a series of ligands | [Interaction fingerprints](fingerprint.md), and [KLIFS kinase pockets](klifs.md) for kinases |
| compare the interactions of several ligands, or feed them to a clustering | [Interaction matrix](matrix.md) |
| use the same positions for different GPCRs (3x32, 6x48...) | [GPCR pockets (GPCRdb)](gpcrdb.md) |
| analyze an entry of the Protein Data Bank without downloading it by hand | [Downloading from the RCSB](rcsb.md) |
| understand why a result differs from another program | [Validation and limits](validation.md) |
| use a file that does not open | [Input formats](formats.md) and [Troubleshooting](troubleshooting.md) |

## Conventions

Commands are shown for a POSIX shell (bash, zsh); on Windows use the same options in a terminal where `ligint.exe` is in the `PATH`.
Distances are in angstroms (A), angles in degrees. A *ligand* is a `HETATM` residue of the file; *protein* is everything else that
is not water. Menu items and buttons of the graphical interface are written in **bold**.
