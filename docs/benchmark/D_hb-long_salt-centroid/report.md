# Benchmark ligint contro PLIP

Data: 2026-10-14. PDB: 50, siti confrontati: 74 (una copia per ligando e PDB).

Sensibilità = quante interazioni di PLIP ritrova ligint. Precisione = quante di ligint sono anche in PLIP. PLIP non è la verità: un disaccordo è un motivo per guardare il caso.

## soglie di ligint

### gruppo D

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 206 | 80 | 45 |  82% |  72% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 32 | 9 | 24 |  57% |  78% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 2 | 2 | 3 |  40% |  50% |
| **totale** | 531 | 232 | 139 |  79% |  70% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 206 | 80 | 45 |  82% |  72% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 32 | 9 | 24 |  57% |  78% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 2 | 2 | 3 |  40% |  50% |
| **totale** | 531 | 232 | 139 |  79% |  70% |

## soglie alla PLIP (ponti salini 5,5 A, alogeni 4,0 A)

### gruppo D

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 206 | 80 | 45 |  82% |  72% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 32 | 9 | 24 |  57% |  78% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 3 | 2 | 2 |  60% |  60% |
| **totale** | 532 | 232 | 138 |  79% |  70% |

### gruppo TUTTI

| tipo | TP | FP | FN | sensibilità | precisione |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 206 | 80 | 45 |  82% |  72% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 32 | 9 | 24 |  57% |  78% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 3 | 2 | 2 |  60% |  60% |
| **totale** | 532 | 232 | 138 |  79% |  70% |

## Cause probabili dei disaccordi (soglie di ligint, tutti i gruppi)

| tipo | causa | n |
|---|---|---:|
| hydrophobic | ligint: contatto idrofobico non riportato da PLIP | 51 |
| hydrophobic | PLIP: contatto idrofobico non riportato da ligint | 35 |
| hbond | ligint: contatto forte (<=3,0 A) non riportato da PLIP | 38 |
| hbond | ligint: contatto 3,0-3,5 A non riportato da PLIP | 34 |
| hbond | PLIP: D-A entro 3,5 A, da esaminare | 24 |
| hbond | PLIP: D-A oltre 3,5 A (fuori soglia di ligint) | 21 |
| hbond | ligint: partner istidina (tautomero assegnato da PLIP) | 8 |
| water | ligint: ponte d'acqua non riportato da PLIP (filtri propri di PLIP) | 82 |
| water | PLIP: ponte d'acqua non riportato da ligint | 29 |
| salt | PLIP: oltre 4,0 A (fuori soglia di ligint) | 18 |
| salt | ligint: ponte salino non riportato da PLIP | 9 |
| salt | PLIP: entro 4,0 A, da esaminare | 6 |
| pistack | ligint: stacking non riportato da PLIP | 3 |
| pistack | PLIP: stacking non riportato da ligint | 1 |
| picat | ligint: cation-pi non riportato da PLIP (protonazione del ligando?) | 5 |
| picat | PLIP: cation-pi non riportato | 2 |
| halogen | PLIP: alogeno non riportato da ligint (F non trattato, soglia) | 3 |
| halogen | ligint: alogeno non riportato | 2 |

## Esito per struttura (soglie di ligint)

| PDB | gruppo | sito | interazioni PLIP | TP | FP | FN |
|---|---|---|---:|---:|---:|---:|
| 1AGM | D | GLC:D:1 | 12 | 8 | 4 | 4 |
| 1AGM | D | MAN:A:485 | 5 | 3 | 2 | 2 |
| 1ATJ | D | HEM:A:350 | 23 | 21 | 4 | 2 |
| 1BMF | D | ADP:D:600 | 18 | 13 | 2 | 5 |
| 1BMF | D | ANP:A:600 | 10 | 9 | 8 | 1 |
| 1BVR | D | NAD:A:1301 | 20 | 16 | 5 | 4 |
| 1BVR | D | THT:A:1302 | 4 | 4 | 5 | 0 |
| 1CKP | D | EDO:A:300 | 1 | 1 | 3 | 0 |
| 1CKP | D | PVB:A:299 | 3 | 3 | 1 | 0 |
| 1EQG | D | BOG:A:801 | 1 | 1 | 1 | 0 |
| 1EQG | D | HEM:A:601 | 16 | 11 | 2 | 5 |
| 1EQG | D | IBP:A:701 | 7 | 7 | 0 | 0 |
| 1EQG | D | NAG:A:1672 | 1 | 1 | 0 | 0 |
| 1FKB | D | RAP:A:108 | 10 | 10 | 2 | 0 |
| 1GIA | D | GSP:A:355 | 25 | 20 | 8 | 5 |
| 1GSS | D | LEE:A:208 | 8 | 7 | 4 | 1 |
| 1HFC | D | PLH:A:280 | 6 | 6 | 3 | 0 |
| 1HPV | D | 478:B:200 | 19 | 17 | 3 | 2 |
| 1I7G | D | AZ2:A:469 | 6 | 5 | 8 | 1 |
| 1I7G | D | CPQ:A:470 | 6 | 6 | 1 | 0 |
| 1ICE | D | ACE:T:285 | 1 | 1 | 1 | 0 |
| 1J91 | D | TBS:A:338 | 1 | 0 | 4 | 1 |
| 1KV2 | D | B96:A:391 | 16 | 12 | 1 | 4 |
| 1KZN | D | CBN:A:1 | 13 | 8 | 5 | 5 |
| 1LDM | D | NAD:A:330 | 14 | 12 | 8 | 2 |
| 1LDM | D | OXM:A:331 | 6 | 5 | 2 | 1 |
| 1LST | D | LYS:A:240 | 12 | 9 | 3 | 3 |
| 1MBN | D | OH:A:154 | 11 | 0 | 1 | 11 |
| 1MCP | D | SO4:H:223 | 3 | 3 | 1 | 0 |
| 1N45 | D | HEM:A:300 | 15 | 12 | 6 | 3 |
| 1N45 | D | SO4:A:234 | 4 | 3 | 1 | 1 |
| 1P0I | D | BUA:A:606 | 6 | 5 | 0 | 1 |
| 1P0I | D | GOL:A:603 | 2 | 2 | 1 | 0 |
| 1P0I | D | MES:A:607 | 5 | 2 | 7 | 3 |
| 1P0I | D | NAG:A:536 | 5 | 1 | 0 | 4 |
| 1P0I | D | SO4:A:601 | 3 | 3 | 0 | 0 |
| 1PTY | D | PTR:A:323 | 15 | 14 | 5 | 1 |
| 1T3R | D | 017:A:1200 | 14 | 12 | 3 | 2 |
| 1T3R | D | PO4:A:502 | 3 | 3 | 1 | 0 |
| 1TPH | D | PGH:1:250 | 11 | 0 | 16 | 11 |
| 1U72 | D | MTX:A:188 | 12 | 10 | 5 | 2 |
| 1U72 | D | NDP:A:187 | 21 | 16 | 7 | 5 |
| 1UDT | D | VIA:A:1000 | 8 | 7 | 2 | 1 |
| 1V4S | D | GLC:A:500 | 11 | 9 | 3 | 2 |
| 1V4S | D | MRK:A:501 | 7 | 7 | 1 | 0 |
| 1XWS | D | BI1:A:1001 | 10 | 9 | 0 | 1 |
| 1XWS | D | SO4:A:2001 | 3 | 3 | 0 | 0 |
| 1YDT | D | IQB:E:351 | 10 | 10 | 5 | 0 |
| 2AQU | D | DR7:B:300 | 18 | 15 | 5 | 3 |
| 2LBD | D | REA:A:500 | 10 | 9 | 3 | 1 |
| 2OHU | D | DMS:A:389 | 3 | 2 | 1 | 1 |
| 2OHU | D | IP7:A:390 | 7 | 7 | 5 | 0 |
| 2VT4 | D | D10:A:406 | 4 | 3 | 0 | 1 |
| 2VT4 | D | P32:A:400 | 9 | 8 | 5 | 1 |
| 2VT4 | D | SOG:A:402 | 5 | 4 | 4 | 1 |
| 3CS9 | D | NIL:A:600 | 16 | 14 | 3 | 2 |
| 3ERD | D | ACY:B:717 | 3 | 3 | 0 | 0 |
| 3ERD | D | DES:A:600 | 13 | 13 | 4 | 0 |
| 3LXK | D | MI1:A:1125 | 6 | 4 | 3 | 2 |
| 3NXU | D | DMS:A:700 | 2 | 1 | 3 | 1 |
| 3NXU | D | HEM:A:508 | 29 | 16 | 3 | 13 |
| 3PBL | D | ETQ:A:1200 | 8 | 6 | 0 | 2 |
| 3PBL | D | GLC:C:1 | 7 | 4 | 7 | 3 |
| 3POZ | D | 03P:A:1023 | 9 | 7 | 8 | 2 |
| 3POZ | D | SO4:A:1 | 2 | 2 | 0 | 0 |
| 3UON | D | BGC:A:1163 | 4 | 4 | 0 | 0 |
| 3UON | D | QNB:A:1162 | 13 | 12 | 3 | 1 |
| 4HG7 | D | NUT:A:201 | 7 | 6 | 3 | 1 |
| 4HG7 | D | SO4:A:202 | 4 | 2 | 3 | 2 |
| 4TMN | D | 0PK:E:317 | 10 | 9 | 3 | 1 |
| 4ZAU | D | YY3:A:1101 | 5 | 5 | 1 | 0 |
| 5L2I | D | LQQ:A:900 | 8 | 6 | 0 | 2 |
| 5P21 | D | GNP:A:167 | 21 | 19 | 9 | 2 |
| 6CEL | D | NAG:A:435 | 4 | 3 | 1 | 1 |

## Problemi

- 1B2Y: nessun ligando a piccola molecola per PLIP
- 1PLC: nessun ligando a piccola molecola per PLIP
- 1URP: nessun ligando a piccola molecola per PLIP
