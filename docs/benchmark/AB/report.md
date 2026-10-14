# Benchmark ligint contro PLIP

Data: 2026-10-14. PDB: 50, siti confrontati: 84 (una copia per ligando e PDB).

Sensibilità = quante interazioni di PLIP ritrova ligint. Precisione = quante di ligint sono anche in PLIP. PLIP non è la verità: un disaccordo è un motivo per guardare il caso.

## soglie di ligint

### gruppo A

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 122 | 21 | 8 |  94% |  85% |
| hbond | 68 | 28 | 25 |  73% |  71% |
| water | 23 | 35 | 14 |  62% |  40% |
| salt | 19 | 1 | 12 |  61% |  95% |
| pistack | 8 | 4 | 1 |  89% |  67% |
| picat | 4 | 3 | 0 | 100% |  57% |
| halogen | 1 | 1 | 2 |  33% |  50% |
| **totale** | 245 | 93 | 62 |  80% |  72% |

### gruppo B

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 161 | 40 | 10 |  94% |  80% |
| hbond | 127 | 20 | 42 |  75% |  86% |
| water | 30 | 43 | 14 |  68% |  41% |
| salt | 17 | 2 | 10 |  63% |  89% |
| pistack | 6 | 3 | 0 | 100% |  67% |
| picat | 1 | 0 | 0 | 100% | 100% |
| halogen | 2 | 0 | 2 |  50% | 100% |
| **totale** | 344 | 108 | 78 |  82% |  76% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 283 | 61 | 18 |  94% |  82% |
| hbond | 195 | 48 | 67 |  74% |  80% |
| water | 53 | 78 | 28 |  65% |  40% |
| salt | 36 | 3 | 22 |  62% |  92% |
| pistack | 14 | 7 | 1 |  93% |  67% |
| picat | 5 | 3 | 0 | 100% |  62% |
| halogen | 3 | 1 | 4 |  43% |  75% |
| **totale** | 589 | 201 | 140 |  81% |  75% |

## soglie alla PLIP (ponti salini 5,5 A, alogeni 4,0 A)

### gruppo A

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 122 | 21 | 8 |  94% |  85% |
| hbond | 68 | 28 | 25 |  73% |  71% |
| water | 23 | 35 | 14 |  62% |  40% |
| salt | 25 | 8 | 6 |  81% |  76% |
| pistack | 8 | 4 | 1 |  89% |  67% |
| picat | 4 | 3 | 0 | 100% |  57% |
| halogen | 1 | 1 | 2 |  33% |  50% |
| **totale** | 251 | 100 | 56 |  82% |  72% |

### gruppo B

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 161 | 40 | 10 |  94% |  80% |
| hbond | 127 | 20 | 42 |  75% |  86% |
| water | 30 | 43 | 14 |  68% |  41% |
| salt | 19 | 8 | 8 |  70% |  70% |
| pistack | 6 | 3 | 0 | 100% |  67% |
| picat | 1 | 0 | 0 | 100% | 100% |
| halogen | 2 | 1 | 2 |  50% |  67% |
| **totale** | 346 | 115 | 76 |  82% |  75% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 283 | 61 | 18 |  94% |  82% |
| hbond | 195 | 48 | 67 |  74% |  80% |
| water | 53 | 78 | 28 |  65% |  40% |
| salt | 44 | 16 | 14 |  76% |  73% |
| pistack | 14 | 7 | 1 |  93% |  67% |
| picat | 5 | 3 | 0 | 100% |  62% |
| halogen | 3 | 2 | 4 |  43% |  60% |
| **totale** | 597 | 215 | 132 |  82% |  74% |

## Cause probabili dei disaccordi (soglie di ligint, tutti i gruppi)

| tipo | causa | n |
|---|---|---:|
| hydrophobic | ligint: contatto idrofobico non riportato da PLIP | 61 |
| hydrophobic | PLIP: contatto idrofobico non riportato da ligint | 18 |
| hbond | PLIP: D-A oltre 3,5 A (fuori soglia di ligint) | 49 |
| hbond | ligint: contatto forte (<=3,0 A) non riportato da PLIP | 26 |
| hbond | PLIP: D-A entro 3,5 A, da esaminare | 18 |
| hbond | ligint: contatto 3,0-3,5 A non riportato da PLIP | 11 |
| hbond | ligint: partner istidina (tautomero assegnato da PLIP) | 11 |
| water | ligint: ponte d'acqua non riportato da PLIP (filtri propri di PLIP) | 78 |
| water | PLIP: ponte d'acqua non riportato da ligint | 28 |
| salt | PLIP: oltre 4,0 A (fuori soglia di ligint) | 17 |
| salt | PLIP: entro 4,0 A, da esaminare | 5 |
| salt | ligint: ponte salino non riportato da PLIP | 3 |
| pistack | ligint: stacking non riportato da PLIP | 7 |
| pistack | PLIP: stacking non riportato da ligint | 1 |
| picat | ligint: cation-pi non riportato da PLIP (protonazione del ligando?) | 3 |
| halogen | PLIP: alogeno non riportato da ligint (F non trattato, soglia) | 4 |
| halogen | ligint: alogeno non riportato | 1 |

## Esito per struttura (soglie di ligint)

| PDB | gruppo | sito | interazioni PLIP | TP | FP | FN |
|---|---|---|---:|---:|---:|---:|
| 181L | B | BNZ:A:400 | 5 | 5 | 0 | 0 |
| 181L | B | HED:A:170 | 0 | 0 | 1 | 0 |
| 1A9U | B | SB2:A:800 | 11 | 10 | 1 | 1 |
| 1AKE | A | AP5:A:215 | 28 | 27 | 12 | 1 |
| 1ANF | B | GLC:B:1 | 13 | 10 | 4 | 3 |
| 1ATP | A | ATP:E:355 | 13 | 11 | 6 | 2 |
| 1AVD | B | BTN:A:400 | 11 | 10 | 2 | 1 |
| 1AVD | B | NAG:A:600 | 2 | 2 | 0 | 0 |
| 1CBS | B | REA:A:200 | 13 | 13 | 0 | 0 |
| 1COY | B | AND:A:508 | 12 | 10 | 2 | 2 |
| 1COY | B | FAD:A:510 | 23 | 18 | 8 | 5 |
| 1CX2 | A | HEM:A:682 | 14 | 9 | 1 | 5 |
| 1CX2 | A | NAG:A:661 | 3 | 0 | 0 | 3 |
| 1CX2 | A | S58:A:701 | 8 | 4 | 2 | 4 |
| 1DB1 | A | VDX:A:428 | 13 | 13 | 5 | 0 |
| 1ERE | A | EST:A:600 | 9 | 9 | 4 | 0 |
| 1FK9 | B | EFZ:A:999 | 7 | 7 | 1 | 0 |
| 1FKJ | A | FK5:A:108 | 8 | 8 | 3 | 0 |
| 1FKN | A | 1OL:C:4 | 6 | 5 | 1 | 1 |
| 1FTJ | B | GLU:A:274 | 13 | 10 | 3 | 3 |
| 1H9Z | A | MYR:A:1001 | 7 | 6 | 3 | 1 |
| 1H9Z | A | RWF:A:3001 | 8 | 6 | 3 | 2 |
| 1HCK | B | ATP:A:400 | 9 | 9 | 4 | 0 |
| 1HEW | A | NAG:B:1 | 10 | 8 | 3 | 2 |
| 1HRC | B | ACE:A:0 | 0 | 0 | 0 | 0 |
| 1HRC | B | HEC:A:105 | 18 | 15 | 8 | 3 |
| 1HVR | B | XK2:A:263 | 18 | 15 | 1 | 3 |
| 1HVY | A | BME:A:1514 | 1 | 0 | 1 | 1 |
| 1HVY | A | D16:A:414 | 12 | 9 | 6 | 3 |
| 1HVY | A | UMP:A:314 | 11 | 10 | 3 | 1 |
| 1HXB | A | ROC:A:100 | 21 | 16 | 4 | 5 |
| 1ICT | A | T44:C:128 | 9 | 6 | 1 | 3 |
| 1M17 | A | AQ4:A:999 | 7 | 6 | 1 | 1 |
| 1M2Z | B | BOG:A:501 | 4 | 2 | 0 | 2 |
| 1M2Z | B | DEX:A:301 | 10 | 7 | 2 | 3 |
| 1QS4 | B | 100:A:1004 | 7 | 4 | 2 | 3 |
| 1R9O | B | FLP:A:501 | 8 | 8 | 2 | 0 |
| 1R9O | B | GOL:A:502 | 2 | 0 | 0 | 2 |
| 1R9O | B | HEM:A:500 | 16 | 14 | 7 | 2 |
| 1RBP | A | RTL:A:183 | 9 | 8 | 2 | 1 |
| 1SA0 | B | CN2:B:700 | 7 | 3 | 0 | 4 |
| 1SA0 | B | GDP:B:602 | 14 | 10 | 0 | 4 |
| 1SA0 | B | GTP:A:600 | 13 | 9 | 1 | 4 |
| 1T46 | A | PO4:A:4 | 4 | 2 | 1 | 2 |
| 1T46 | A | STI:A:3 | 19 | 16 | 5 | 3 |
| 1U19 | B | ACE:A:0 | 0 | 0 | 0 | 0 |
| 1U19 | B | HTG:A:1507 | 2 | 2 | 2 | 0 |
| 1U19 | B | HTO:B:1401 | 3 | 2 | 1 | 1 |
| 1U19 | B | NAG:C:1 | 2 | 2 | 0 | 0 |
| 1U19 | B | PLM:A:1322 | 5 | 5 | 2 | 0 |
| 1U19 | B | RET:A:1296 | 7 | 7 | 2 | 0 |
| 1US0 | B | CIT:A:400 | 8 | 4 | 3 | 4 |
| 1US0 | B | LDT:A:320 | 8 | 6 | 6 | 2 |
| 1US0 | B | NDP:A:318 | 27 | 25 | 2 | 2 |
| 1UW6 | A | NCT:A:1208 | 4 | 3 | 9 | 1 |
| 1UWH | A | BAX:A:1723 | 12 | 8 | 2 | 4 |
| 1VRT | A | NVP:A:999 | 8 | 8 | 3 | 0 |
| 1YET | B | GDM:A:275 | 13 | 11 | 5 | 2 |
| 2BAT | A | NAG:A:470 | 3 | 2 | 0 | 1 |
| 2BAT | A | SIA:A:600 | 13 | 10 | 2 | 3 |
| 2CPP | B | CAM:A:422 | 3 | 3 | 2 | 0 |
| 2CPP | B | HEM:A:417 | 16 | 11 | 4 | 5 |
| 2HHB | A | HEM:A:142 | 13 | 12 | 3 | 1 |
| 2HHB | A | PO4:B:147 | 0 | 0 | 0 | 0 |
| 2ITY | B | IRE:A:2020 | 4 | 3 | 4 | 1 |
| 2PRG | A | BRL:A:1 | 8 | 6 | 3 | 2 |
| 2RH1 | B | 12P:A:416 | 1 | 0 | 0 | 1 |
| 2RH1 | B | ACM:A:411 | 2 | 1 | 1 | 1 |
| 2RH1 | B | BU1:A:409 | 7 | 6 | 1 | 1 |
| 2RH1 | B | CAU:A:408 | 10 | 10 | 3 | 0 |
| 2RH1 | B | CLR:A:412 | 7 | 7 | 1 | 0 |
| 2RH1 | B | GLC:B:1 | 5 | 3 | 0 | 2 |
| 2RH1 | B | PLM:A:415 | 0 | 0 | 0 | 0 |
| 2RH1 | B | SO4:A:402 | 5 | 4 | 0 | 1 |
| 3CPA | A | GLY:A:501 | 13 | 6 | 1 | 7 |
| 3EML | A | SO4:A:407 | 3 | 2 | 1 | 1 |
| 3EML | A | STE:A:402 | 2 | 2 | 1 | 0 |
| 3EML | A | ZMA:A:401 | 8 | 7 | 1 | 1 |
| 3ERT | B | OHT:A:600 | 13 | 11 | 4 | 2 |
| 3MXF | B | DMS:A:171 | 5 | 2 | 0 | 3 |
| 3MXF | B | EDO:A:170 | 4 | 4 | 0 | 0 |
| 3MXF | B | JQ1:A:1 | 10 | 9 | 2 | 1 |
| 4DFR | B | MTX:A:161 | 12 | 8 | 7 | 4 |
| 5P9J | B | 8E8:A:701 | 7 | 7 | 7 | 0 |

## Problemi

- 3CLN: nessun ligando a piccola molecola per PLIP
- 8ADH: nessun ligando a piccola molecola per PLIP
