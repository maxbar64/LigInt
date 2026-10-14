# Benchmark of ligint against PLIP

Date: 2026-10-14. PDB files: 50, sites compared: 74 (one copy per ligand and PDB).

Sensitivity = how many PLIP interactions ligint finds. Precision = how many of ligint's are also in PLIP. PLIP is not the truth: a disagreement is a reason to look at the case.

## ligint thresholds

### group D

| type | TP | FP | FN | sensitivity | precision |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 181 | 61 | 70 |  72% |  75% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 29 | 7 | 27 |  52% |  81% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 2 | 2 | 3 |  40% |  50% |
| **total** | 503 | 211 | 167 |  75% |  70% |

### group ALL

| type | TP | FP | FN | sensitivity | precision |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 181 | 61 | 70 |  72% |  75% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 29 | 7 | 27 |  52% |  81% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 2 | 2 | 3 |  40% |  50% |
| **total** | 503 | 211 | 167 |  75% |  70% |

## PLIP-like thresholds (salt bridges 5.5 A, halogens 4.0 A)

### group D

| type | TP | FP | FN | sensitivity | precision |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 181 | 61 | 70 |  72% |  75% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 33 | 20 | 23 |  59% |  62% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 3 | 2 | 2 |  60% |  60% |
| **total** | 508 | 224 | 162 |  76% |  69% |

### group ALL

| type | TP | FP | FN | sensitivity | precision |
|---|---:|---:|---:|---:|---:|
| hydrophobic | 217 | 51 | 35 |  86% |  81% |
| hbond | 181 | 61 | 70 |  72% |  75% |
| water | 55 | 82 | 29 |  65% |  40% |
| salt | 33 | 20 | 23 |  59% |  62% |
| pistack | 15 | 3 | 1 |  94% |  83% |
| picat | 4 | 5 | 2 |  67% |  44% |
| halogen | 3 | 2 | 2 |  60% |  60% |
| **total** | 508 | 224 | 162 |  76% |  69% |

## Probable causes of the disagreements (ligint thresholds, all groups)

| type | cause | n |
|---|---|---:|
| hydrophobic | ligint: hydrophobic contact not reported by PLIP | 51 |
| hydrophobic | PLIP: hydrophobic contact not reported by ligint | 35 |
| hbond | PLIP: D-A beyond 3.5 A (outside the ligint threshold) | 44 |
| hbond | ligint: strong contact (<=3.0 A) not reported by PLIP | 39 |
| hbond | PLIP: D-A within 3.5 A, to be examined | 26 |
| hbond | ligint: 3.0-3.5 A contact not reported by PLIP | 16 |
| hbond | ligint: histidine partner (tautomer assigned by PLIP) | 6 |
| water | ligint: water bridge not reported by PLIP (PLIP applies its own filters) | 82 |
| water | PLIP: water bridge not reported by ligint | 29 |
| salt | PLIP: beyond 4.0 A (outside the ligint threshold) | 21 |
| salt | ligint: salt bridge not reported by PLIP | 7 |
| salt | PLIP: within 4.0 A, to be examined | 6 |
| pistack | ligint: stacking not reported by PLIP | 3 |
| pistack | PLIP: stacking not reported by ligint | 1 |
| picat | ligint: cation-pi not reported by PLIP (ligand protonation?) | 5 |
| picat | PLIP: cation-pi not reported by ligint | 2 |
| halogen | PLIP: halogen bond not reported by ligint (F not handled, threshold) | 3 |
| halogen | ligint: halogen bond not reported by PLIP | 2 |

## Outcome per structure (ligint thresholds)

| PDB | group | site | PLIP interactions | TP | FP | FN |
|---|---|---|---:|---:|---:|---:|
| 1AGM | D | GLC:D:1 | 12 | 7 | 3 | 5 |
| 1AGM | D | MAN:A:485 | 5 | 3 | 2 | 2 |
| 1ATJ | D | HEM:A:350 | 23 | 20 | 4 | 3 |
| 1BMF | D | ADP:D:600 | 18 | 13 | 2 | 5 |
| 1BMF | D | ANP:A:600 | 10 | 9 | 7 | 1 |
| 1BVR | D | NAD:A:1301 | 20 | 14 | 5 | 6 |
| 1BVR | D | THT:A:1302 | 4 | 4 | 4 | 0 |
| 1CKP | D | EDO:A:300 | 1 | 1 | 2 | 0 |
| 1CKP | D | PVB:A:299 | 3 | 3 | 1 | 0 |
| 1EQG | D | BOG:A:801 | 1 | 1 | 1 | 0 |
| 1EQG | D | HEM:A:601 | 16 | 10 | 2 | 6 |
| 1EQG | D | IBP:A:701 | 7 | 7 | 0 | 0 |
| 1EQG | D | NAG:A:1672 | 1 | 1 | 0 | 0 |
| 1FKB | D | RAP:A:108 | 10 | 10 | 2 | 0 |
| 1GIA | D | GSP:A:355 | 25 | 18 | 8 | 7 |
| 1GSS | D | LEE:A:208 | 8 | 7 | 4 | 1 |
| 1HFC | D | PLH:A:280 | 6 | 5 | 3 | 1 |
| 1HPV | D | 478:B:200 | 19 | 13 | 3 | 6 |
| 1I7G | D | AZ2:A:469 | 6 | 5 | 8 | 1 |
| 1I7G | D | CPQ:A:470 | 6 | 6 | 1 | 0 |
| 1ICE | D | ACE:T:285 | 1 | 1 | 1 | 0 |
| 1J91 | D | TBS:A:338 | 1 | 0 | 4 | 1 |
| 1KV2 | D | B96:A:391 | 16 | 12 | 1 | 4 |
| 1KZN | D | CBN:A:1 | 13 | 8 | 5 | 5 |
| 1LDM | D | NAD:A:330 | 14 | 10 | 6 | 4 |
| 1LDM | D | OXM:A:331 | 6 | 4 | 2 | 2 |
| 1LST | D | LYS:A:240 | 12 | 9 | 3 | 3 |
| 1MBN | D | OH:A:154 | 11 | 0 | 1 | 11 |
| 1MCP | D | SO4:H:223 | 3 | 3 | 1 | 0 |
| 1N45 | D | HEM:A:300 | 15 | 11 | 6 | 4 |
| 1N45 | D | SO4:A:234 | 4 | 2 | 1 | 2 |
| 1P0I | D | BUA:A:606 | 6 | 5 | 0 | 1 |
| 1P0I | D | GOL:A:603 | 2 | 2 | 1 | 0 |
| 1P0I | D | MES:A:607 | 5 | 1 | 7 | 4 |
| 1P0I | D | NAG:A:536 | 5 | 1 | 0 | 4 |
| 1P0I | D | SO4:A:601 | 3 | 3 | 0 | 0 |
| 1PTY | D | PTR:A:323 | 15 | 14 | 5 | 1 |
| 1T3R | D | 017:A:1200 | 14 | 12 | 3 | 2 |
| 1T3R | D | PO4:A:502 | 3 | 3 | 1 | 0 |
| 1TPH | D | PGH:1:250 | 11 | 0 | 13 | 11 |
| 1U72 | D | MTX:A:188 | 12 | 9 | 5 | 3 |
| 1U72 | D | NDP:A:187 | 21 | 16 | 6 | 5 |
| 1UDT | D | VIA:A:1000 | 8 | 7 | 2 | 1 |
| 1V4S | D | GLC:A:500 | 11 | 8 | 3 | 3 |
| 1V4S | D | MRK:A:501 | 7 | 7 | 0 | 0 |
| 1XWS | D | BI1:A:1001 | 10 | 8 | 0 | 2 |
| 1XWS | D | SO4:A:2001 | 3 | 3 | 0 | 0 |
| 1YDT | D | IQB:E:351 | 10 | 10 | 5 | 0 |
| 2AQU | D | DR7:B:300 | 18 | 15 | 5 | 3 |
| 2LBD | D | REA:A:500 | 10 | 9 | 3 | 1 |
| 2OHU | D | DMS:A:389 | 3 | 1 | 1 | 2 |
| 2OHU | D | IP7:A:390 | 7 | 6 | 5 | 1 |
| 2VT4 | D | D10:A:406 | 4 | 3 | 0 | 1 |
| 2VT4 | D | P32:A:400 | 9 | 7 | 3 | 2 |
| 2VT4 | D | SOG:A:402 | 5 | 3 | 3 | 2 |
| 3CS9 | D | NIL:A:600 | 16 | 14 | 3 | 2 |
| 3ERD | D | ACY:B:717 | 3 | 3 | 0 | 0 |
| 3ERD | D | DES:A:600 | 13 | 13 | 3 | 0 |
| 3LXK | D | MI1:A:1125 | 6 | 4 | 2 | 2 |
| 3NXU | D | DMS:A:700 | 2 | 1 | 3 | 1 |
| 3NXU | D | HEM:A:508 | 29 | 15 | 3 | 14 |
| 3PBL | D | ETQ:A:1200 | 8 | 6 | 0 | 2 |
| 3PBL | D | GLC:C:1 | 7 | 3 | 4 | 4 |
| 3POZ | D | 03P:A:1023 | 9 | 7 | 7 | 2 |
| 3POZ | D | SO4:A:1 | 2 | 2 | 0 | 0 |
| 3UON | D | BGC:A:1163 | 4 | 3 | 0 | 1 |
| 3UON | D | QNB:A:1162 | 13 | 12 | 3 | 1 |
| 4HG7 | D | NUT:A:201 | 7 | 6 | 3 | 1 |
| 4HG7 | D | SO4:A:202 | 4 | 2 | 3 | 2 |
| 4TMN | D | 0PK:E:317 | 10 | 9 | 3 | 1 |
| 4ZAU | D | YY3:A:1101 | 5 | 5 | 1 | 0 |
| 5L2I | D | LQQ:A:900 | 8 | 6 | 0 | 2 |
| 5P21 | D | GNP:A:167 | 21 | 19 | 8 | 2 |
| 6CEL | D | NAG:A:435 | 4 | 3 | 1 | 1 |

## Problems

- 1B2Y: no small-molecule ligand for PLIP
- 1PLC: no small-molecule ligand for PLIP
- 1URP: no small-molecule ligand for PLIP
