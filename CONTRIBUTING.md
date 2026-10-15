# Contributing

## Build and test

```
cmake -S . -B build -G Ninja -DLIGINT_WARNINGS_AS_ERRORS=ON
cmake --build build
QT_QPA_PLATFORM=offscreen ctest --test-dir build --output-on-failure
```

Before changing the library it is also worth checking the sanitizers:
`cmake --preset sanitize && cmake --build --preset sanitize && ctest --preset sanitize`.

## Project rules

- **All text in the code is in English** (comments, messages, GUI labels, documentation, test messages).
- Every bug fix starts from a **test that reproduces it** (`tests/run_tests.py` for the chemical rules and the
  parsing, `tests/unit/test_core.cpp` for the API). The synthetic tests use fragments with known geometry: see the
  molecule builders in `run_tests.py`.
- New thresholds go into `ligint::Params`, into `Params::set` (same name as the command-line option), into the
  `--help` text and into the README. If the GUI exposes the setting, add it to `MainWindow::currentArgs` (the command line
  shown in the Log): the GUI smoke test fails when the generated command is not equivalent to the window.
- A new fingerprint class needs an entry in `FpClass`, `fpClassName`/`fpClassTitle`, a glyph in `src/core/svg.cpp`, the color in
  `gui/Style.cpp` (`fpClassColor`) and the documentation in `docs/FINGERPRINT.md`.
- The interaction types are the constants in `ligint::type` (also used by the GUI in `gui/Style.cpp` and by the
  comparison scripts): when adding one, update the library, `gui/Style.cpp` and `tests/compare_plip.py`.
- The library does not print: warnings and information go into `Result`.
- Changes that alter the results must be measured with `tests/benchmark.py` on structures that were **not** used to
  decide the change. A comparison with PLIP measures the agreement with PLIP, not correctness.

## Network features

The RCSB download and the KLIFS service are optional (libcurl, json-c). The tests of group `net` never use the network: transfers that must fail are sent to a
proxy that does not exist (`https_proxy`), the others use the cache or a local mapping file. Keep it so: a test that reaches a real server is not reproducible. Build once
with and once without the libraries, and keep the exported CMake package free of build-only targets (the CI builds a program against the installed package).

## The manual

The manual is built from Markdown (see `docs/ARCHITECTURE.md`, section *The manual*). After editing `docs/manual/*.md` or one of the documents that
it includes, run `pip install -r tools/requirements-manual.txt` once and `python3 tools/build_manual.py`, and commit the changes in `docs/manual/html/`: the
CI fails when they are out of date. A new page goes into `docs/manual/manifest.json`. Interface labels, option names and threshold defaults quoted in the
manual are checked against the program by the tests; keep them in sync.

## Installers and releases

The packages are described in `cmake/Packaging.cmake` and `cmake/CPackOptions.cmake` (CPack) and tried by the `linux-packages`
and `windows` jobs of the CI. See `docs/INSTALL.md` for how to build them and how to publish a release (tag `v<version>`).
A new installed file needs a component (`COMPONENT Runtime` or `Development`) or it ends up in no package.

## Benchmark

```
pip install plip openbabel-wheel
python3 tests/benchmark.py folder_with_pdb/ --groups groups.txt --out results/
```

The script runs PLIP (reusing the reports already in `results/plip/`) and ligint, and writes `report.md`,
`disagreements.csv` and `per_category.csv`. The groups used so far are in `docs/benchmark/groups.tsv`.
