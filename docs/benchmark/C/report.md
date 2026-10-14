# Benchmark ligint contro PLIP

Data: 2026-10-14. PDB: 50, siti confrontati: 83 (una copia per ligando e PDB).

Sensibilità = quante interazioni di PLIP ritrova ligint. Precisione = quante di ligint sono anche in PLIP. PLIP non è la verità: un disaccordo è un motivo per guardare il caso.

## soglie di ligint

### gruppo C

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 242 | 60 | 15 |  94% |  80% |
| hbond | 176 | 54 | 70 |  72% |  77% |
| water | 66 | 97 | 35 |  65% |  40% |
| salt | 33 | 11 | 32 |  51% |  75% |
| pistack | 21 | 8 | 1 |  95% |  72% |
| picat | 5 | 5 | 1 |  83% |  50% |
| halogen | 3 | 0 | 8 |  27% | 100% |
| **totale** | 546 | 235 | 162 |  77% |  70% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 242 | 60 | 15 |  94% |  80% |
| hbond | 176 | 54 | 70 |  72% |  77% |
| water | 66 | 97 | 35 |  65% |  40% |
| salt | 33 | 11 | 32 |  51% |  75% |
| pistack | 21 | 8 | 1 |  95% |  72% |
| picat | 5 | 5 | 1 |  83% |  50% |
| halogen | 3 | 0 | 8 |  27% | 100% |
| **totale** | 546 | 235 | 162 |  77% |  70% |

## soglie alla PLIP (ponti salini 5,5 A, alogeni 4,0 A)

### gruppo C

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 242 | 60 | 15 |  94% |  80% |
| hbond | 176 | 54 | 70 |  72% |  77% |
| water | 66 | 97 | 35 |  65% |  40% |
| salt | 39 | 22 | 26 |  60% |  64% |
| pistack | 21 | 8 | 1 |  95% |  72% |
| picat | 5 | 5 | 1 |  83% |  50% |
| halogen | 6 | 1 | 5 |  55% |  86% |
| **totale** | 555 | 247 | 153 |  78% |  69% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 242 | 60 | 15 |  94% |  80% |
| hbond | 176 | 54 | 70 |  72% |  77% |
| water | 66 | 97 | 35 |  65% |  40% |
| salt | 39 | 22 | 26 |  60% |  64% |
| pistack | 21 | 8 | 1 |  95% |  72% |
| picat | 5 | 5 | 1 |  83% |  50% |
| halogen | 6 | 1 | 5 |  55% |  86% |
| **totale** | 555 | 247 | 153 |  78% |  69% |

## Cause probabili dei disaccordi (soglie di ligint, tutti i gruppi)

| tipo | causa | n |
|---|---|---:|
| hydrophobic | ligint: contatto idrofobico non riportato da PLIP | 60 |
| hydrophobic | PLIP: contatto idrofobico non riportato da ligint | 15 |
| hbond | PLIP: D-A oltre 3,5 A (fuori soglia di ligint) | 50 |
| hbond | ligint: contatto forte (<=3,0 A) non riportato da PLIP | 27 |
| hbond | PLIP: D-A entro 3,5 A, da esaminare | 20 |
| hbond | ligint: contatto 3,0-3,5 A non riportato da PLIP | 17 |
| hbond | ligint: partner istidina (tautomero assegnato da PLIP) | 10 |
| water | ligint: ponte d'acqua non riportato da PLIP (filtri propri di PLIP) | 97 |
| water | PLIP: ponte d'acqua non riportato da ligint | 35 |
| salt | PLIP: oltre 4,0 A (fuori soglia di ligint) | 29 |
| salt | ligint: ponte salino non riportato da PLIP | 11 |
| salt | PLIP: entro 4,0 A, da esaminare | 3 |
| pistack | ligint: stacking non riportato da PLIP | 8 |
| pistack | PLIP: stacking non riportato da ligint | 1 |
| picat | ligint: cation-pi non riportato da PLIP (protonazione del ligando?) | 5 |
| picat | PLIP: cation-pi non riportato | 1 |
| halogen | PLIP: alogeno non riportato da ligint (F non trattato, soglia) | 8 |

## Esito per struttura (soglie di ligint)

| PDB | gruppo | sito | interazioni PLIP | TP | FP | FN |
|---|---|---|---:|---:|---:|---:|
| 1A28 | C | STR:A:1 | 10 | 9 | 3 | 1 |
| 1A3K | C | NAG:B:1 | 9 | 7 | 3 | 2 |
| 1A6M | C | HEM:A:154 | 19 | 15 | 4 | 4 |
| 1A6M | C | SO4:A:155 | 3 | 2 | 2 | 1 |
| 1ABE | C | ARA:A:307 | 8 | 4 | 5 | 4 |
| 1ACJ | C | THA:A:999 | 11 | 10 | 4 | 1 |
| 1AJ6 | C | NOV:A:1 | 9 | 7 | 3 | 2 |
| 1AQ1 | C | STU:A:299 | 9 | 8 | 2 | 1 |
| 1ATK | C | E64:A:216 | 12 | 9 | 4 | 3 |
| 1BSX | C | T3:A:1 | 7 | 4 | 5 | 3 |
| 1C3W | C | LI1:A:601 | 5 | 5 | 0 | 0 |
| 1C3W | C | RET:A:301 | 9 | 9 | 3 | 0 |
| 1C3W | C | SQU:A:701 | 3 | 3 | 1 | 0 |
| 1D3G | C | ACT:A:401 | 4 | 1 | 0 | 3 |
| 1D3G | C | BRE:A:397 | 14 | 12 | 2 | 2 |
| 1D3G | C | DDQ:A:700 | 3 | 3 | 1 | 0 |
| 1D3G | C | FMN:A:398 | 18 | 15 | 3 | 3 |
| 1D3G | C | ORO:A:399 | 9 | 8 | 2 | 1 |
| 1D3G | C | SO4:A:400 | 3 | 2 | 2 | 1 |
| 1DWC | C | MIT:H:1 | 14 | 9 | 3 | 5 |
| 1E3G | C | R18:A:1000 | 11 | 10 | 5 | 1 |
| 1E8X | C | ATP:A:3000 | 10 | 4 | 3 | 6 |
| 1FJS | C | GOL:A:502 | 4 | 3 | 0 | 1 |
| 1FJS | C | Z34:A:500 | 11 | 8 | 7 | 3 |
| 1GD1 | C | NAD:O:336 | 19 | 14 | 7 | 5 |
| 1GD1 | C | SO4:O:338 | 4 | 4 | 1 | 0 |
| 1HWK | C | 117:A:2 | 14 | 13 | 3 | 1 |
| 1HWK | C | ADP:A:101 | 4 | 3 | 1 | 1 |
| 1HXW | C | RIT:B:301 | 17 | 15 | 3 | 2 |
| 1KIF | C | BEZ:A:349 | 5 | 5 | 1 | 0 |
| 1KIF | C | FAD:A:348 | 24 | 19 | 6 | 5 |
| 1O86 | C | GLY:A:2000 | 2 | 1 | 0 | 1 |
| 1O86 | C | LPR:A:702 | 17 | 9 | 7 | 8 |
| 1OYA | C | FMN:A:401 | 19 | 16 | 3 | 3 |
| 1PAU | C | ACE:C:501 | 1 | 1 | 0 | 0 |
| 1PAU | C | ASJ:C:505 | 8 | 4 | 1 | 4 |
| 1PPB | C | 0G6:H:1 | 11 | 8 | 5 | 3 |
| 1QCF | C | PP1:A:532 | 11 | 10 | 1 | 1 |
| 1QSG | C | GLC:A:1303 | 1 | 1 | 1 | 0 |
| 1QSG | C | NAD:A:1301 | 19 | 16 | 4 | 3 |
| 1QSG | C | TCL:A:1302 | 7 | 6 | 2 | 1 |
| 1S9J | C | ATP:A:535 | 8 | 6 | 6 | 2 |
| 1S9J | C | BBM:A:1001 | 7 | 5 | 4 | 2 |
| 1T69 | C | SHH:A:379 | 6 | 3 | 4 | 3 |
| 1VID | C | DNC:A:302 | 3 | 3 | 4 | 0 |
| 1VID | C | SAM:A:301 | 16 | 10 | 7 | 6 |
| 1X70 | C | 715:A:801 | 11 | 6 | 3 | 5 |
| 1X70 | C | NDG:C:2 | 0 | 0 | 4 | 0 |
| 1XKK | C | FMM:A:91 | 13 | 12 | 4 | 1 |
| 1XKK | C | PO4:A:81 | 2 | 1 | 0 | 1 |
| 2A65 | C | BOG:A:701 | 6 | 5 | 1 | 1 |
| 2A65 | C | LEU:A:601 | 7 | 5 | 5 | 2 |
| 2CYP | C | HEM:A:295 | 20 | 17 | 1 | 3 |
| 2FDN | C | SF4:A:61 | 1 | 0 | 0 | 1 |
| 2GQG | C | 1N1:A:501 | 10 | 10 | 3 | 0 |
| 2GQG | C | GOL:B:1 | 2 | 1 | 2 | 1 |
| 2HU4 | C | G39:A:800 | 8 | 6 | 2 | 2 |
| 2IFB | C | PLM:A:133 | 9 | 7 | 3 | 2 |
| 2VCI | C | 2GJ:A:1225 | 11 | 10 | 7 | 1 |
| 2W26 | C | RIV:A:1001 | 4 | 4 | 3 | 0 |
| 2XP2 | C | VGH:A:9000 | 6 | 5 | 0 | 1 |
| 2ZD1 | C | EDO:B:429 | 2 | 1 | 3 | 1 |
| 2ZD1 | C | SO4:A:556 | 3 | 2 | 0 | 1 |
| 2ZD1 | C | T27:A:557 | 12 | 11 | 2 | 1 |
| 3GRS | C | FAD:A:479 | 23 | 15 | 14 | 8 |
| 3GRS | C | PO4:A:480 | 6 | 3 | 4 | 3 |
| 3KS3 | C | GOL:A:304 | 2 | 2 | 1 | 0 |
| 3LN1 | C | BOG:A:703 | 1 | 1 | 2 | 0 |
| 3LN1 | C | CEL:A:682 | 14 | 10 | 4 | 4 |
| 3LN1 | C | HEM:A:605 | 15 | 12 | 2 | 3 |
| 3LON | C | BME:A:998 | 0 | 0 | 2 | 0 |
| 3LON | C | NNA:A:999 | 14 | 11 | 5 | 3 |
| 3MEC | C | 65B:A:561 | 11 | 11 | 4 | 0 |
| 3MEC | C | SO4:A:562 | 3 | 2 | 0 | 1 |
| 3OG7 | C | 032:A:1 | 14 | 12 | 5 | 2 |
| 4AZU | C | NO3:A:900 | 2 | 2 | 0 | 0 |
| 4DKL | C | 1PE:A:617 | 0 | 0 | 0 | 0 |
| 4DKL | C | BF0:A:601 | 11 | 8 | 4 | 3 |
| 4DKL | C | CLR:A:614 | 7 | 7 | 1 | 0 |
| 4DKL | C | MPG:A:615 | 6 | 2 | 0 | 4 |
| 4DKL | C | SO4:A:602 | 5 | 2 | 2 | 3 |
| 4G5J | C | 0WM:A:1102 | 1 | 0 | 3 | 1 |
| 4G5J | C | 0WN:A:1101 | 8 | 4 | 1 | 4 |

## Problemi

- 2HNP: nessun ligando a piccola molecola per PLIP
