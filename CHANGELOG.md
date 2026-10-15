# Changelog

## 0.5.1

Installers and packages (no change to the programs):

- CPack configuration: `ligint` and `ligint-dev` `.deb` packages and a portable `.tar.gz` for Linux (`.rpm` with
  `cpack -G RPM`), an NSIS installer and a portable `.zip` for Windows. The installers include the programs, the
  documentation, the manual pages (`man ligint`, `man ligint-gui`), the menu entry with icon and the PDB file
  association (Linux), the Qt libraries (Windows), the Start menu entry and an optional `PATH` update (Windows).
- Application icon, window icon and Windows version information.
- CI: `linux-packages` and `windows` jobs build, try (install, run, uninstall) and keep the installers as artifacts;
  a `v*` tag publishes them in a GitHub release.
- Install rules split in two components (`Runtime`, `Development`); `docs/INSTALL.md`.
- The Windows installer and the MSVC build have not been run by the authors: the first CI run will show whether they need adjustments.

## 0.5.0

- `--altloc X` (library: `Options::altloc`) chooses the alternate conformation: for every residue with alternate locations the
  conformation X is used if the residue has it, otherwise its first one (as before). The program reports for how many
  residues X was used and warns if X does not exist in the first model. Occupancies are not considered.
- The listing of the ligands now also shows the alternate conformations present in the file
  (`Alternate conformations present: A, B`; `Result::altlocs`); it appears only for files that have them, so the listing of
  those files changes (27 of the 119 structures used in the tests), while the rest of the output is unchanged
  (1075 runs compared with 0.4.0).
- GUI: *Alternate conformation* selector in the Options group (filled from the file), included in the copied command line.
- Tests: C++ unit tests, 2 new Python tests (also on a ligand with alternate conformations, ENDMDL, invalid arguments, speed on
  200,000 atoms) and the GUI smoke test (combo, command line and equivalence with the CLI for the chosen conformation).

## 0.4.0

- GUI: the equivalent `ligint` command line of each analysis is written in the Log tab and can be copied
  (*File > Copy command line*, *Copy fingerprint command line*). Only non-default settings appear; the PLIP-like preset is
  written as the explicit options it sets; quoting is for POSIX shells.
- The command line selects a ligand by residue number with an optional insertion code: `501` (any), `501A` (code A),
  `501-` (no code); a malformed residue number is now an error instead of being read as 0. In the library: `Options::icode`.
- The ligand header and `LigandResult` now show/carry the insertion code (`=== Ligand LIG 5A chain L ...`); the output is
  otherwise unchanged (1075 runs compared with 0.3.0).
- Tests: the GUI smoke test runs the generated commands through `/bin/sh` with the real `ligint` and compares them with the
  result shown by the window (default settings, preset and options, SDF template, insertion codes, fingerprint SVG, shell quoting).

## 0.3.0

Interaction fingerprints (one row per protein residue, one column per interaction class):

- new library API in `ligint/fingerprint.hpp`: `makeFingerprint`, `PocketDefinition` (plain-text file mapping residues
  to a position label and a region), CSV and JSON writers, and `fingerprintSvg` (figure with the residues grouped by
  region and one glyph per class, legend and tooltips; no Qt needed);
- the interaction records now say which side donates the hydrogen and the charge of the residue
  (`HitRecord::hbResidueRole`, `HitRecord::residueCharge`), which the fingerprint classes need;
- command line: `--fingerprint csv|json|svg` and `--pocket FILE`;
- GUI: **Fingerprint** tab (figure and table), pocket definition loading, export to SVG/CSV/JSON; the figure needs the
  optional Qt Svg module, otherwise the GUI shows the table only; the control panel scrolls on small windows;
- documentation in `docs/FINGERPRINT.md` with an illustrative pocket definition (`docs/examples/`);
- 6 new Python tests, new C++ unit tests and GUI smoke test checks for fingerprints and pocket files.

The existing output of `ligint` is unchanged.

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
