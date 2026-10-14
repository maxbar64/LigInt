# Changelog

## 0.2.0

All text in the code is now in English: comments, messages, GUI labels, documentation. **This changes the program
output and is not compatible with 0.1.0:**

- interaction types are now English (`Hydrogen bond`, `Salt bridge`, `Hydrophobic contact`, `pi-stacking (parallel)`,
  `pi-stacking (T-shaped)`, `Cation-pi`, `Halogen bond`, `Metal coordination`, `Water bridge`, `Steric clash`, ...) and
  are available as constants in `ligint::type`;
- the CSV header is `ligand,type,favorable,ligand_atom,protein_atom,distance,notes`;
- table headings, summaries, warnings and error messages are English; ring labels use `ring` instead of `anello`;
- rows are sorted by type name, so their order differs from 0.1.0 (the content is identical: 28,861 CSV rows over 956
  runs were compared, after mapping the texts, with no differences);
- the benchmark scripts write English reports and CSV columns (`group`, `site`, `type`, ...).

Fixed:

- the coordination sphere of a metal no longer omits the protein nitrogens bound to it (His in zinc sites, introduced by
  the 0.1.0 fix that made them non-acceptors), with a regression test;
- the code now compiles with Clang (a lambda captured a structured binding, which is a C++20 extension);
- the CMake package exports the library as `ligint::core`.

## 0.1.0

First version organized as a CMake project.

- `ligint_core` library with a public API (`analyzeFile`, `analyzeStream`, `printText`, `printCsv`, `loadSdf`, `loadMol2`);
  the result includes the coordinates of the interactions, the ligand atoms and the atoms of the residues in contact
  for the graphics.
- Command-line program `ligint`: output identical to that of the single-file version (verified on 1075 runs of 119
  complexes in several modes), plus `--version`.
- Qt GUI `ligint-gui` with 3D view, tables, filters by type, CSV and image export, SDF/MOL2 templates.
- Interactions: hydrogen bonds (also long ones with a possible geometry), salt bridges (also by centroids),
  hydrophobic, pi-stacking, cation-pi, halogens, metals, water bridges, clashes; multi-residue ligands merged;
  DNA/RNA as partners.
- Tests: C++ (`core_unit`), Python, GUI smoke test (offscreen), comparison and benchmark against PLIP.
