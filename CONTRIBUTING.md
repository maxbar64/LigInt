# Contribuire

## Compilare e provare

```
cmake -S . -B build -G Ninja -DLIGINT_WARNINGS_AS_ERRORS=ON
cmake --build build
QT_QPA_PLATFORM=offscreen ctest --test-dir build --output-on-failure
```

Prima di una modifica alla libreria conviene controllare anche i sanitizer:
`cmake --preset sanitize && cmake --build --preset sanitize && ctest --preset sanitize`.

## Regole del progetto

- Ogni correzione di un difetto parte da un **test che lo riproduce** (`tests/run_tests.py` per le regole chimiche e
  il parsing, `tests/unit/test_core.cpp` per l'API). I test sintetici usano frammenti con geometria nota: vedi i
  costruttori di molecole in `run_tests.py`.
- Le soglie nuove vanno aggiunte a `ligint::Params`, a `Params::set` (nome uguale all'opzione da riga di comando),
  all'`--help` e al README.
- I tipi di interazione sono stringhe in italiano usate anche dalla GUI (`gui/Style.cpp`) e dagli script di confronto:
  aggiungendone uno vanno aggiornati tutti e tre.
- La libreria non stampa: avvisi e informazioni vanno in `Result`.
- Le modifiche che cambiano i risultati vanno misurate con `tests/benchmark.py` su strutture **non** usate per
  decidere la modifica. Un confronto con PLIP misura l'accordo con PLIP, non la correttezza.

## Benchmark

```
pip install plip openbabel-wheel
python3 tests/benchmark.py cartella_con_pdb/ --groups gruppi.txt --out risultati/
```

Lo script esegue PLIP (riusando i report già presenti in `risultati/plip/`) e ligint, e scrive `report.md`,
`disagreements.csv` e `per_category.csv`. I gruppi usati finora sono in `docs/benchmark/groups.tsv`.
