# Architettura

## Livelli

```
include/ligint/ligint.hpp   API pubblica: Params, Options, Result, loadSdf/loadMol2, analyzeFile/analyzeStream, printText/printCsv
src/core/
  internal.hpp              tipi interni (Atom, Role, Hit, Ring, Grid) e funzioni condivise tra i file
  analysis.cpp              Params, lettura del PDB e analisi (analyzeStream)
  chemistry.cpp             anelli e ruoli chimici (donatore, accettore, carica) di proteina e ligando
  templates.cpp             lettura SDF/MOL2 e abbinamento al ligando per coordinate
  report.cpp                rapporti testuali (tabella e CSV)
src/cli/main.cpp            riga di comando: argomenti -> Options -> analyzeFile -> printText/printCsv
gui/                        Qt: MainWindow, ViewerWidget (vista 3D), HitsModel/HitsFilterProxy, Style
```

La libreria non dipende da Qt e non stampa: tutto ciò che il programma comunica (avvisi, informazioni, risultati)
è in `Result`. La GUI e la CLI sono due client della stessa API.

## Flusso di `analyzeStream`

1. **Lettura del PDB** (primo modello): atomi pesanti, conformazione alternativa (la prima incontrata per residuo),
   coordinate non valide scartate, acque a parte, H del PDB assegnati all'atomo pesante vicino (se ≥ 10 per la proteina).
   Se non è richiesto alcun ligando (o non esiste) il risultato è solo l'elenco (`listingOnly`).
2. **Ruoli della proteina:** donatore/accettore/carica/idrofobico per residuo e atomo; raffinati con gli H espliciti;
   N/O coordinati a un metallo non fanno H-bond.
3. **Ligandi composti:** i residui HETATM legati covalentemente (non metalli) sono uniti in un solo ligando.
4. **Per ogni ligando:** connettività da template (SDF/MOL2, per coordinate) o geometrica; ruoli del ligando
   (`ligandRole`); coordinazione ai metalli; ricerca atomo-atomo su una griglia spaziale (H-bond, ponti salini,
   idrofobici, alogeni, metalli, clash); ponti salini per centroidi (opzionali); anelli (π-stacking, cation-π);
   ponti d'acqua.
5. **Risultato:** interazioni ordinate, con coordinate (per la grafica), atomi e legami del ligando, atomi dei residui a
   contatto e sfera di coordinazione dei metalli.

## Aggiungere un tipo di interazione

1. Calcolarlo in `analysis.cpp` (nel ciclo atomo-atomo o dopo, come per anelli e ponti d'acqua) con una nuova stringa di tipo
   e le coordinate `p1`/`p2` (e `pw` per l'acqua).
2. Dargli uno stile in `gui/Style.cpp` e, se va confrontato con PLIP, mapparlo in `tests/compare_plip.py`.
3. Scrivere il test sintetico in `tests/run_tests.py` e documentare la soglia (`Params`, `Params::set`, `--help`, README).

## GUI

`MainWindow` possiede un `HitsModel` (tabella) con `HitsFilterProxy` (filtro per tipo, ordinamento numerico), un `ViewerWidget` e
le tabelle dei residui e del registro. L'analisi gira con `QtConcurrent::run` e `QFutureWatcher`; `analyze(true)` la esegue in modo
sincrono (usato dai test). `ViewerWidget` fa il rendering in software con `QPainter` (proiezione ortografica, ordinamento per
profondità, orientazione iniziale sugli assi principali) e può disegnare su un `QImage` (`renderImage`), usato per salvare le
immagini e nei test.

## Test

| Livello | Dove | Cosa |
|---|---|---|
| API | `tests/unit/test_core.cpp` | soglie, rapporti, template, ligandi composti, input malformati |
| Chimica e parsing | `tests/run_tests.py` | frammenti sintetici con esito noto, SDF/MOL2, altloc, modelli, input anomali, 9 complessi reali in `tests/data/` |
| GUI | `tests/gui/gui_smoke.cpp` | apertura, analisi, filtri, rendering e salvataggio (piattaforma offscreen) |
| Accordo con PLIP | `tests/compare_plip.py`, `tests/benchmark.py` | misura contro uno strumento esterno (non in CTest) |
