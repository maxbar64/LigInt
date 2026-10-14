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
  `--help` text and into the README.
- The interaction types are the constants in `ligint::type` (also used by the GUI in `gui/Style.cpp` and by the
  comparison scripts): when adding one, update the library, `gui/Style.cpp` and `tests/compare_plip.py`.
- The library does not print: warnings and information go into `Result`.
- Changes that alter the results must be measured with `tests/benchmark.py` on structures that were **not** used to
  decide the change. A comparison with PLIP measures the agreement with PLIP, not correctness.

## Benchmark

```
pip install plip openbabel-wheel
python3 tests/benchmark.py folder_with_pdb/ --groups groups.txt --out results/
```

The script runs PLIP (reusing the reports already in `results/plip/`) and ligint, and writes `report.md`,
`disagreements.csv` and `per_category.csv`. The groups used so far are in `docs/benchmark/groups.tsv`.
