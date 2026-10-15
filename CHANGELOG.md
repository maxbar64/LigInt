# Changelog

## 0.8.2

- **KLIFS: the pocket is a MOL2 file.** The second real run showed that the request for the pocket answers with a `pocket.mol2` (the residues of the pocket, not a table), which the
  program could not read. It now reads it: the residues are the positions 1 to 85 that the structure has, the regions are those of the KLIFS alignment, and the 85-character pocket
  sequence (`_` for a missing residue; the field `pocket` of the service, or `--klifs-sequence`) places them when some are missing. The sequence is checked against the file
  (number of residues, and every standard residue against its letter), so a wrong one is refused. `--klifs-mapping` accepts a `pocket.mol2` too (offline, GUI included).
  Checked on four real KLIFS pocket files (BMX, BRAF, ABL1 chains A and B): the missing residues agree with KLIFS, and the ABL1 hinge falls on Glu316-Phe317-Met318.
- The colors of the regions did not apply to the names that KLIFS uses: `αC`, `αD`, `αE`, `αF` (Greek letters) and `xDFG`.
- Errors of the online path say what the answer contained when the sequence is missing; the structure list is now requested also when `--klifs-id` is given, to get the sequence.
- Tests: the reading of the MOL2 in C++ and from the command line (all the error cases, corruption, the Greek region names through the SVG figure).
- Still not confirmed against the live service: the name of the sequence field in the answer about the structures and the alternate-model field; a complete online run on a structure with missing residues.

## 0.8.1

- **KLIFS service, first real run.** It reached the service and found two structures for the same chain (IDs 1445 and 1453), but the alternate model was read as empty,
  so the program could not choose and stopped. Now: the alternate model is read from `alt`, `alternate_model`, `structure.alternate_model`, `alt_model` or `model`;
  the error says to choose with `--klifs-id N` (or `--klifs-model M` when the answer names the model), and when the model is not in the answer it lists the fields that
  the answer has; a response that is not JSON, an empty one, or a pocket row without a position are reported with the beginning of the answer or the names of the fields.
- **GUI:** when several KLIFS structures match, a window lists the candidates (ID, chain, model) and the request is repeated with the one that you pick.
  Library: `KlifsAmbiguousError` (with the candidates) and `selectKlifsStructure`, tested with answers in the format of the service.
- Still not confirmed: the name of the alternate-model field and the answer of the pocket request (`structure_get_pocket`).

## 0.8.0

- **Downloads from the RCSB** (`--download-pdb ID`, `--download-format`, `--download-output`, `--download-cache`, `--download-refresh`, `--download-only`;
  GUI: *File > Download from RCSB...*) and **KLIFS kinase pocket definitions** (`--klifs-pdb`, `--klifs-chain`, `--klifs-mapping`, `--klifs-output`,
  `--klifs-model`, `--klifs-id`; GUI: *Generate from KLIFS...*), in C++ with libcurl and json-c, both optional at build time. New library headers
  `ligint/rcsb.hpp` and `ligint/klifs.hpp`; `tools/klifs_to_pocket.py` is a Python alternative for the online step. Manual pages: *Downloading from the RCSB*,
  *KLIFS kinase pockets*.
- **Fixed in review:** with network support the installed CMake package could not be used (`find_package(ligint)` failed on a build-only target; it now names
  the system libraries, and the `ligint-dev` and `ligint-devel` packages depend on their development packages); `--download-pdb` could not be combined with
  a ligand name (the arguments after the options are now the ligand); `--klifs-*` silently lost the pocket when there was no structure and no `--klifs-output`,
  and ignored it without `--fingerprint` (an error and a warning); the interface blocked for up to 90 s during a download or a KLIFS request (they run now
  in a background thread, like the analysis) and offered the online KLIFS source in builds that cannot use it; libcurl older than 7.85 (RHEL 9) did not
  compile; the KLIFS requests did not restrict the protocol of redirects.
- Tests: group `net` (KLIFS mapping variants, errors and corrupted files; the pocket in the fingerprint; the command-line rules; the RCSB cache, the ID
  validation and the failed transfers, with no real network traffic), the GUI background tasks, a C++ test of the RCSB cache.
- CI: the Linux jobs install libcurl and json-c, so that the network code is built and tested; the Rocky 9 job tries them too (libcurl 7.76).
- Not verified: a real transfer from the RCSB and from the KLIFS service (the test environment cannot reach them); the address `https://klifs.net/api_v2` used by
  the native KLIFS client is unconfirmed.

## 0.7.0

- **User manual in HTML, inside the program.** *Help > User manual* (`F1`) opens a window with the manual: table of contents,
  full-text search over all the pages, back/forward history, zoom, find in the page, links between pages. `F1` opens the page of the tab
  you are looking at (Interactions, Residues, Fingerprint, Log); *Help > Command-line reference* goes straight to the options. The pages are
  embedded in the executable (no files to install, works offline) and are also installed as plain HTML in `share/doc/ligint/manual/`
  (*Manual in web browser*, or open `index.html`).
- Fourteen pages: welcome, quick start (with a real output), the graphical interface, the command-line reference (every option and
  threshold), interactions and assumptions, input formats, fingerprints, validation and limits, comparison with PLIP, installation,
  troubleshooting, developers, changelog, about. Numbers that appear in several documents come from a single source (the manual includes
  sections of the README).
- The manual is generated from Markdown by `tools/build_manual.py` (the HTML is committed, so building ligint still needs no Python).
  The CI checks that the HTML is up to date.
- Tests: the committed HTML is up to date with its sources; every link, anchor and image of the pages resolves; the thresholds,
  defaults and options in the command-line reference are the ones of `ligint --help`; the example output of the quick start is the real output of the
  program; the GUI test opens the manual, navigates, searches and checks that every context-help target (`F1`) exists.
- Man pages and README point to the manual.

## 0.6.0

- **PDBx/mmCIF input.** `ligint structure.cif LIG` gives the same results as the PDB file of the same entry: the format is recognized
  from the content (`--format pdb|cif|auto`, library: `Options::format`, `Result::format`). The `_atom_site` category is read as a stream
  (first model, author chain and numbering, `label_*` as fallback), so very large structures are fine. mmCIF allows what the PDB
  format cannot: residue names of five characters (new ligands) and chains of several characters.
- Residue labels with a chain of several characters are written `ASP 189 [AA]` (the labels of one-character chains do not change);
  the fingerprint uses the same rule, also for the positions of a pocket definition without contacts.
- gzip-compressed files are refused with a clear message instead of being read as garbage.
- The readers are separated from the analysis (`pdb.cpp`, `mmcif.cpp` produce the same internal atom record): the output for PDB files
  is unchanged (1314 runs on 119 structures, also with `--altloc`, compared with 0.5.1).
- GUI: opens `.cif` and `.mmcif` files (dialog, window texts).
- Tests: 119 real structures converted to mmCIF with gemmi give the same results as the PDB files (1314 runs, 0 differences);
  C++ and Python tests for quoting, text fields, comments, null values, long names and chains, models, alternate conformations,
  malformed and randomly corrupted files, speed on 200,000 atoms; `tests/data/1AZM.cif` (gemmi conversion) in the GUI test.
- `docs/INPUT_FORMATS.md`: columns used, labels, limits. The reader has not been tried on files downloaded from the PDB archive.

## 0.5.1

Installers and packages (no change to the programs):

- CPack configuration: `ligint` and `ligint-dev` `.deb` packages and a portable `.tar.gz` for Linux (`.rpm` with
  `cpack -G RPM`), an NSIS installer and a portable `.zip` for Windows. The installers include the programs, the
  documentation, the manual pages (`man ligint`, `man ligint-gui`), the menu entry with icon and the PDB file
  association (Linux), the Qt libraries (Windows), the Start menu entry and an optional `PATH` update (Windows).
- Application icon, window icon and Windows version information.
- CI: `linux-packages` and `windows` jobs build, try (install, run, uninstall) and keep the installers as artifacts;
  `rpm-rocky9` builds the RPM in a Rocky Linux 9 container and `rpm-alma9` installs it (and the development package) in a clean
  AlmaLinux 9 container; a `v*` tag publishes everything in a GitHub release. The workflow can also be started by hand.
- Install rules split in two components (`Runtime`, `Development`); `docs/INSTALL.md`.
- The Windows installer, the MSVC build and the RPM jobs have not been run by the authors: the first CI run will show whether they need adjustments.

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
