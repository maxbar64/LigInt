# Architecture

## Layers

```
include/ligint/ligint.hpp   public API: Params, Options, Result, loadSdf/loadMol2, analyzeFile/analyzeStream, printText/printCsv, ligint::type
src/core/
  internal.hpp              internal types (Atom, Role, Hit, Ring, Grid) and functions shared between the files
  analysis.cpp              Params and the analysis (analyzeStream): applies the rules to the atoms that the readers produce
  pdb.cpp, mmcif.cpp        readers of the PDB and PDBx/mmCIF formats: both produce RawAtom (internal.hpp), nothing else differs
  chemistry.cpp             rings and chemical roles (donor, acceptor, charge) of protein and ligand
  templates.cpp             SDF/MOL2 reading and matching to the ligand by coordinates
  report.cpp                text reports (table and CSV)
  fingerprint.cpp           pocket definition, per-residue classification, CSV/JSON fingerprint output (ligint/fingerprint.hpp)
  svg.cpp                   SVG figure of a fingerprint
src/cli/main.cpp            command line: arguments -> Options -> analyzeFile -> printText/printCsv
gui/                        Qt: MainWindow, ViewerWidget (3D view), SvgView (figure), HelpWindow (manual), HitsModel/HitsFilterProxy, Style
```

The library does not depend on Qt and does not print: everything the program communicates (warnings, information,
results) is in `Result`. The GUI and the CLI are two clients of the same API.

## Flow of `analyzeStream`

1. **Reading** (PDB or mmCIF, recognized from the content; first model): heavy atoms, alternative conformation (the first encountered per residue), invalid
   coordinates discarded, waters kept apart, PDB hydrogens assigned to the nearby heavy atom (if there are at least 10 for the protein).
   If no ligand is requested (or it does not exist) the result is only the listing (`listingOnly`).
2. **Protein roles:** donor/acceptor/charge/hydrophobic per residue and atom; refined with explicit hydrogens;
   N/O coordinated to a metal make no H-bonds.
3. **Composite ligands:** covalently linked HETATM residues (not metals) are merged into a single ligand.
4. **For each ligand:** connectivity from a template (SDF/MOL2, by coordinates) or geometric; ligand roles
   (`ligandRole`); coordination to metals; atom-atom search on a spatial grid (H-bonds, salt bridges,
   hydrophobic, halogens, metals, clashes); salt bridges by centroids (optional); rings (π-stacking, cation-π);
   water bridges.
5. **Result:** sorted interactions, with coordinates (for the graphics), ligand atoms and bonds, atoms of the residues in
   contact and the coordination sphere of the metals.

## The manual

The user manual is Markdown: the files of `docs/manual/` and some of the other documents (`docs/FINGERPRINT.md`, `docs/INPUT_FORMATS.md`,
`docs/INSTALL.md`, `docs/BENCHMARK.md`, `docs/ARCHITECTURE.md`, `CHANGELOG.md`), listed in `docs/manual/manifest.json`. A page can include a section of
another file with `<!-- include: README.md "Section title" -->`, so that a number (a threshold, a table) has a single source.
`python3 tools/build_manual.py` converts everything to HTML in `docs/manual/html/` (pages, images, `toc.json`, `manual.qrc`); that
folder is **committed**, so building ligint needs no Python. `HelpWindow` (`gui/HelpWindow.cpp`) shows the pages from the Qt resources
of `manual.qrc` in a `QTextBrowser` (contents tree from `toc.json`, full-text search over the plain text of the pages); the files are also installed in
`share/doc/ligint/manual/`. The markup is kept simple because the Qt text browser understands only a subset of HTML and CSS: tables get `border`
attributes, headings get `<a name>` anchors, images get a `width`. The CI and the tests check that the HTML is up to date and that its links resolve.

## Fingerprints

`makeFingerprint` works only on a `LigandResult`: it classifies each favorable `HitRecord` (type, `hbResidueRole`,
`residueCharge`) into the `FpClass` columns and groups the counts by residue; a `PocketDefinition` fixes which residues
are listed and in which order. The SVG and the GUI table are views of the same `Fingerprint`. See `docs/FINGERPRINT.md`.

## Adding an interaction type

1. Add a constant in `ligint::type` (`ligint.hpp`) and compute the interaction in `analysis.cpp` (inside the atom-atom
   loop or afterwards, as for rings and water bridges) with the coordinates `p1`/`p2` (and `pw` for water).
2. Give it a style in `gui/Style.cpp`, decide its fingerprint class in `classesOf` (`src/core/fingerprint.cpp`) and, if it must be
   compared with PLIP, map it in `tests/compare_plip.py`.
3. Write the synthetic test in `tests/run_tests.py` and document the threshold (`Params`, `Params::set`, `--help`, README).

## GUI

`MainWindow` owns a `HitsModel` (table) with a `HitsFilterProxy` (filter by type, numeric sorting), a `ViewerWidget` and
the tables of residues and log. The analysis runs with `QtConcurrent::run` and `QFutureWatcher`; `analyze(true)` runs it
synchronously (used by the tests). `ViewerWidget` renders in software with `QPainter` (orthographic projection, depth
sorting, initial orientation along the principal axes) and can draw on a `QImage` (`renderImage`), used to save the
images and in the tests.

## Tests

| Level | Where | What |
|---|---|---|
| API | `tests/unit/test_core.cpp` | thresholds, reports, templates, composite ligands, malformed input |
| Chemistry and parsing | `tests/run_tests.py` | synthetic fragments with known outcome, SDF/MOL2, altlocs, models, anomalous input, 9 real complexes in `tests/data/` |
| GUI | `tests/gui/gui_smoke.cpp` | opening, analysis, filters, rendering and saving (offscreen platform) |
| Agreement with PLIP | `tests/compare_plip.py`, `tests/benchmark.py` | measurement against an external tool (not in CTest) |
