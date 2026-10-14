# ligint

Interazioni non covalenti **ligando-proteina** da un file PDB: legami a idrogeno, ponti salini, contatti
idrofobici, π-stacking, cation-π, legami ad alogeno, coordinazione di metalli e ponti d'acqua, più clash e
repulsioni. Il progetto contiene:

| Componente | Cosa è |
|---|---|
| `ligint_core` | libreria C++17 senza dipendenze (API in `include/ligint/ligint.hpp`) |
| `ligint` | programma a riga di comando (tabella leggibile o CSV) |
| `ligint-gui` | interfaccia grafica Qt con vista 3D delle interazioni, tabelle e filtri |
| `tests/` | test C++, suite Python, test di fumo della GUI, confronto con PLIP |

![GUI su 1HSG (indinavir)](docs/img/gui-1HSG-window.png)

## Compilare

Requisiti: **CMake ≥ 3.16**, un compilatore **C++17**; per la GUI **Qt 6** (o Qt 5) con Widgets e Concurrent.
Se Qt non c'è, CMake avvisa e compila solo libreria e riga di comando.

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure      # con QT_QPA_PLATFORM=offscreen se non c'è uno schermo
cmake --install build --prefix /percorso/di/installazione
```

Con CMake ≥ 3.21 si possono usare i preset: `cmake --preset release && cmake --build --preset release`
(altri preset: `debug`, `cli-only`, `sanitize`).

| Opzione CMake | Default | Significato |
|---|---|---|
| `LIGINT_BUILD_CLI` | ON | eseguibile `ligint` |
| `LIGINT_BUILD_GUI` | ON | `ligint-gui` (richiede Qt; se manca, avviso e prosegue) |
| `LIGINT_BUILD_TESTS` | ON | test registrati in CTest |
| `LIGINT_SANITIZE` | OFF | AddressSanitizer + UBSan (GCC/Clang; meglio senza GUI) |
| `LIGINT_WARNINGS_AS_ERRORS` | OFF | avvisi del compilatore come errori |
| `LIGINT_BUILD_SCREENSHOT_TOOL` | OFF | `ligint-screenshot`, per le immagini della documentazione |

## Uso

### Riga di comando

```
ligint struttura.pdb                       elenca i ligandi (HETATM) presenti, con il numero di copie
ligint struttura.pdb ATP                   analizza tutti i ligandi ATP
ligint struttura.pdb ATP A 501             solo catena A, residuo 501
ligint struttura.pdb ATP --csv             output CSV (le informazioni di servizio vanno su stderr)
ligint struttura.pdb LIG --sdf lig.sdf     annotazione chimica da SDF (anche --mol2)
ligint struttura.pdb LIG --plip-like       preset per il confronto con PLIP
ligint --help                              soglie e valori di default
```

Opzioni: `--csv`, `--no-water`, `--no-merge`, `--plip-like`, `--sdf F`, `--mol2 F`, `--metal-donors`,
`--backbone-hyd`, `--hb-long X`, `--salt-centroid X`, `--version` e tutte le soglie (`--hb-max`, `--salt-max`,
`--hyd-max`, …; vedi `--help`). Opzioni sconosciute, soglie non numeriche o file mancanti terminano con errore
(codice 1); soglie oltre 6 Å per le ricerche a griglia vengono limitate con un avviso.

### Interfaccia grafica

```
ligint-gui [struttura.pdb]
```

Apri un PDB (menu, trascinamento o argomento), scegli il ligando, regola opzioni e soglie, premi **Analizza**.
La vista 3D (trascina: ruota; tasto destro o centrale: sposta; rotellina: zoom; doppio clic: ripristina) mostra
ligando, residui a contatto, ioni e acque, con le interazioni colorate per tipo, distanze, etichette e legenda.
Selezionando una riga della tabella l'interazione si evidenzia nella vista. Si possono filtrare i tipi, aggiungere
template SDF/MOL2, esportare le interazioni in CSV e salvare l'immagine in PNG. L'analisi gira in un thread separato.

![Vista 3D su 4IKF (inibitore dell'integrasi, Mg, DNA)](docs/img/gui-4IKF-viewer.png)

### Come libreria

```cpp
#include <ligint/ligint.hpp>

ligint::Options opt;
opt.ligName = "ATP";
opt.params.applyPlipLike();                    // facoltativo
ligint::Result r = ligint::analyzeFile("1ATP.pdb", opt);
for (const auto& lig : r.ligands)
    for (const auto& h : lig.hits)             // tipo, atomi, distanza, coordinate per la grafica...
        std::cout << h.type << " " << h.ligandAtom << " " << h.proteinAtom << " " << h.distance << "\n";
```

La libreria non stampa nulla: avvisi e informazioni sono nel `Result`. `printText` e `printCsv` riproducono
l'output della riga di comando. In un altro progetto CMake, dopo `cmake --install`:

```cmake
find_package(ligint REQUIRED)
target_link_libraries(mio_programma PRIVATE ligint::core)
```

## Cosa riconosce (soglie di default)

| Tipo | Criterio |
|---|---|
| Legame a idrogeno | D···A 2,2–3,5 Å, donatore e accettore compatibili; con H noti D–H···A ≥ 120° e H···A ≤ 2,7 Å, altrimenti angolo antecedente ≥ 90° (ammine a tre sostituenti: asse N–H entro 50°) |
| Ponte salino | cariche opposte, atomi a ≤ 4,0 Å; con H-bond sulla stessa coppia riga unica "anche legame a idrogeno" |
| Repulsione elettrostatica | cariche uguali a ≤ 4,5 Å |
| Contatto idrofobico | C–C apolari a 3,3–4,0 Å (carboni non legati a N/O/P) |
| π-stacking | centroidi ≤ 5,5 Å; parallelo (normali ≤ 30°, piani a 3,0–4,5 Å, offset ≤ 3,5 Å) o a T (≥ 60°, offset ≤ 2,0 Å) |
| Cation-π | catione–centroide ≤ 6,0 Å, entro 30° dall'asse dell'anello |
| Legame ad alogeno | Cl/Br/I···N/O/S ≤ 3,5 Å, angolo C–X···A ≥ 140° (il fluoro non è considerato) |
| Coordinazione metallica | O/N con coppia solitaria o S entro il limite per elemento (Mg 2,5; Zn 2,6; Ca 2,8; …) |
| Ponte d'acqua | acqua a 2,2–3,5 Å da ligando e proteina, angolo ligando–O–proteina 75–140°, senza contatto diretto |
| Clash | distanza < somma dei raggi di van der Waals − 0,5 Å |
| Possibile legame covalente | distanza < 2,0 Å |

### Assunzioni chimiche

I PDB cristallografici non hanno idrogeni, quindi la chimica del ligando è dedotta dalla geometria:

- **pH fisiologico:** carbossilati, fosfati e solfonati (due O terminali sullo stesso atomo) sono anionici; ammine
  alifatiche con soli sostituenti sp3, amidine e guanidini non in anello sono cationiche.
- **Ammidi:** un N legato a un C carbonilico non è accettore. N di anello con coppia solitaria è accettore.
- **Metalli:** O/N coordinati a un metallo sono considerati deprotonati (non donatori) salvo `--metal-donors`;
  His e altri N legati a un metallo non fanno H-bond.
- **Istidina:** donatore e accettore, mai carica.
- **Proteina:** DNA/RNA sono partner (fosfati, basi, stacking); se il PDB contiene gli H della proteina vengono usati.
- **Si legge solo il primo modello;** di ogni residuo si tiene la prima conformazione alternativa incontrata.
- **Ligandi di più residui** (disaccaridi, glicani, dipeptidi) legati covalentemente sono un unico ligando
  (come in PLIP); `--no-merge` li tiene separati.

### SDF e MOL2

`--sdf` / `--mol2` forniscono ordini di legame, idrogeni espliciti e cariche formali. Il ligando del PDB è abbinato al
template **per coordinate** (tolleranza 0,6 Å, stesso numero di atomi pesanti): il file deve avere le stesse
coordinate del PDB. Con idrogeni nel template, donatori/accettori e angoli D–H···A sono reali; senza H, o se
l'abbinamento fallisce, si ricade sull'euristica geometrica (il programma lo dichiara). SDF V3000 non supportato.

### Preset e opzioni per il confronto con PLIP

- `--hb-long 4.1`: H-bond con D···A fino a 4,1 Å solo se l'H (reale o ideale) può puntare all'accettore.
- `--salt-centroid 5.5`: ponti salini tra i centroidi dei gruppi carichi (un contatto per gruppo).
- `--plip-like` = `--hb-long 4.1 --salt-centroid 5.5 --hal-max 4.0`; le opzioni esplicite prevalgono.

Sono spente di default: `--hb-long` è stata confermata su dati non usati per la taratura, `--salt-centroid` ha un
effetto piccolo e incerto.

## Verifiche

```
ctest --test-dir build --output-on-failure           core_unit + python_suite + gui_smoke
python3 tests/run_tests.py [--bin build/src/ligint] [--sanitize]     66 test Python (regole chimiche, SDF/MOL2,
                                                     parsing, input malformati, complessi reali in tests/data)
python3 tests/compare_plip.py *.pdb                  ligint contro PLIP (pip install plip openbabel-wheel)
python3 tests/benchmark.py cartella/ --groups g.txt  pipeline completa con report e CSV dei disaccordi
```

**Accordo con PLIP** (una copia per ligando e PDB; sensibilità = interazioni di PLIP ritrovate; precisione =
interazioni di ligint presenti anche in PLIP; intervalli al 95% con bootstrap sulle strutture):

| Gruppo | Siti | Sensibilità | Precisione |
|---|---:|---|---|
| A, B (usati per le correzioni) | 84 | 81% [77–84] | 75% [71–78] |
| C (cieco, poi esaminato) | 83 | 77% [73–80] | 70% [66–73] |
| D (cieco), default | 47 strutture | 75% [70–80] | 70% [66–74] |
| D, `--plip-like` | 47 strutture | 79% [74–84] | 70% [65–74] |

Dettagli, metodo e limiti in [`docs/BENCHMARK.md`](docs/BENCHMARK.md); report e CSV dei disaccordi in `docs/benchmark/`.
**PLIP non è la verità:** molti disaccordi sono convenzioni (vedi sotto), e il confronto non sostituisce una
validazione con dati sperimentali.

## Limiti noti e differenze da PLIP

- PLIP conta come legami ad alogeno anche quelli del **fluoro**; ligint no.
- PLIP tratta l'**istidina come carica positiva** (circa un terzo dei ponti salini mancanti); ligint no.
- PLIP conta contatti idrofobici di carboni di residui polari e dell'eme; ligint li limita ai residui idrofobici.
- **Ponti d'acqua:** PLIP applica filtri propri (un donatore già in un H-bond non forma un ponte, al massimo due ponti
  da donatore per acqua); ligint no e riporta molti più ponti (precisione ≈ 40%).
- PLIP non analizza acidi nucleici come partner; ligint sì (non validati contro un riferimento).
- Non coperti: SDF V3000, strutture NMR con più modelli, pH diversi da quello fisiologico, legami covalenti con la
  proteina (segnalati solo come "possibile legame covalente"), MOL2/SDF reali da docking (provati solo l'SDF di 1EFY
  e un MOL2 sintetico).

## Stato dei test

Verificato: Linux, GCC 13 e Clang 18 con `-Wall -Wextra -Werror`, Qt 6.4 con piattaforma offscreen (test di fumo e
screenshot), build con ASan+UBSan, installazione e uso come pacchetto CMake. **Non verificato:** macOS e Windows, Qt 5 (il
codice evita API solo-Qt6, ma non è stato compilato). I workflow in `.github/workflows/` **non sono mai stati
eseguiti**: i comandi sono quelli provati qui in locale su Linux (compresi `-Werror`, clang e sanitizer), ma la prima
esecuzione su GitHub potrebbe richiedere ritocchi. Le GUI sono provate in modalità offscreen: l'uso interattivo con mouse e finestre di dialogo non è stato
provato a mano.

## Struttura

```
include/ligint/ligint.hpp   API pubblica
src/core/                   libreria: analysis.cpp, chemistry.cpp, templates.cpp, report.cpp, internal.hpp
src/cli/main.cpp            riga di comando
gui/                        GUI Qt (MainWindow, ViewerWidget, HitsModel, Style) e strumento per gli screenshot
tests/                      unit/ (C++), gui/ (test di fumo), run_tests.py, compare_plip.py, benchmark.py, data/ (PDB di prova)
docs/                       ARCHITECTURE.md, BENCHMARK.md, benchmark/ (report), img/
cmake/                      file di configurazione del pacchetto CMake
```

Altri documenti: [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md), [`CONTRIBUTING.md`](CONTRIBUTING.md), [`CHANGELOG.md`](CHANGELOG.md).

## Licenza

MIT, vedi [`LICENSE`](LICENSE) (il titolare del copyright è un segnaposto: modificalo prima di pubblicare). I file PDB in
`tests/data/` sono dati del Protein Data Bank (dominio pubblico, CC0); `1efy.sdf` è stato generato con Protoss
(ProteinsPlus) a partire dal PDB 1EFY.
