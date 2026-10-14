# Confronto con PLIP

Questo documento descrive come è stato misurato l'accordo tra `ligint` e [PLIP](https://github.com/pharmai/plip)
(versione 3.0.1; sul report 2.4.0 di 1AZM fornito dall'utente dà le stesse interazioni), che cosa si può e non si
può concludere, e come ripetere la misura.

**PLIP non è la verità.** Il confronto misura l'accordo con un altro strumento, che ha convenzioni proprie
(protonazione con Open Babel, soglie diverse, filtri sui ponti d'acqua). Un disaccordo è un motivo per guardare il
caso, non un errore automatico. Non c'è un confronto con dati sperimentali indipendenti.

## Metodo

- **Strutture:** 150 complessi PDB in quattro gruppi (A e B da 25, C e D da 50; elenco in `groups.tsv`), scelti per coprire famiglie
  proteiche e chimiche diverse (chinasi, proteasi, recettori nucleari, GPCR, enzimi con cofattori ed eme, metalloproteine,
  macrocicli, zuccheri, ligandi alogenati...). Gli ID sono stati scelti a memoria e controllati sul titolo dopo il
  download.
- **Siti:** per ogni PDB PLIP analizza tutti i ligandi a piccola molecola; ligint è eseguito sullo stesso ligando
  (chain, residuo). Si conta **una sola copia per ligando e PDB** (altrimenti strutture con 20 copie dominano i totali).
  I ligandi classificati da PLIP come polimeri (catene lunghe di glicani), peptidi e acidi nucleici sono esclusi.
- **Confronto a livello di residuo e tipo:** (tipo, residuo, numero, catena). Tipi confrontati: idrofobico, H-bond, ponte
  d'acqua, ponte salino, π-stacking, cation-π, alogeno. Non confrontati: clash, repulsioni, legami covalenti, contatti con la
  catena principale, coordinazione metallica (PLIP la descrive tra proteina e metallo, ligint tra ligando e metallo).
- **Convenzione sugli H-bond con cariche opposte:** PLIP riporta insieme ponte salino e H-bond sulla stessa coppia solo
  nel 13% dei casi; gli H-bond di ligint "assistiti da carica" (o fusi nel ponte salino) non entrano nel punteggio
  dell'H-bond.
- **Sensibilità** = interazioni di PLIP ritrovate da ligint; **precisione** = interazioni di ligint presenti anche in
  PLIP. Gli intervalli al 95% sono bootstrap sulle strutture (le interazioni della stessa struttura non sono
  indipendenti).

## Gruppi e uso nello sviluppo

| Gruppo | Ruolo |
|---|---|
| A, B | usati per cercare errori e correggere il programma (eme, ponti d'acqua, residui legati, ammine, amidine...) |
| C | primo gruppo non usato per le correzioni; ne sono stati però esaminati i disaccordi |
| D | non usato per nessuna decisione tranne la verifica delle opzioni `--hb-long` e `--salt-centroid` (proposte dopo aver visto A, B e C) |

I risultati di A, B e C sono quindi ottimisti rispetto a complessi mai visti; la stima più onesta è D.

## Risultati (soglie di default, una copia per ligando)

| Gruppo | Sensibilità | Precisione |
|---|---|---|
| A | 80% [73–87] | 72% [68–77] |
| B | 82% [78–85] | 76% [72–81] |
| C | 77% [73–80] | 70% [66–73] |
| D | 75% [70–80] | 70% [66–74] |
| D con `--hb-long 4.1 --salt-centroid 5.5` | 79% [74–84] | 70% [65–74] |

Gruppo D per tipo (default):

| Tipo | Sensibilità | Precisione |
|---|---:|---:|
| idrofobico | 86% | 81% |
| legame a idrogeno | 72% | 75% |
| ponte d'acqua | 65% | 40% |
| ponte salino | 52% | 81% |
| π-stacking | 94% | 83% |
| cation-π | 67% | 44% |
| alogeno | 40% | 50% |

Effetto delle due opzioni su D (differenza appaiata, bootstrap per struttura, opzioni − default):

| | Sensibilità | Precisione |
|---|---|---|
| H-bond (`--hb-long 4.1`) | +9,9 punti [+5,9; +14,8] | −2,7 [−5,0; −0,3] |
| Ponti salini (`--salt-centroid 5.5`) | +5,1 [0,0; +11,8] | −2,3 [−9,3; +2,5] |
| Totale | +4,2 [+2,7; +5,9] | −0,8 [−1,8; +0,1] |

Le opzioni sono quindi spente di default: `--hb-long` ha un effetto netto e confermato, `--salt-centroid` incerto.

## Perché ligint e PLIP divergono

- **Fluoro:** PLIP conta come legami ad alogeno anche quelli con F (12 dei 32 casi dei primi due benchmark); ligint solo Cl, Br, I.
- **Istidina:** PLIP la considera carica positiva; circa il 36% dei ponti salini che mancano a ligint la coinvolge.
- **H-bond lunghi:** circa i due terzi dei contatti che trova solo PLIP stanno tra 3,5 e 4,1 Å (limite di PLIP, con geometria dell'H
  calcolata e vincolo angolare molto più blando); `--hb-long` ne recupera una parte.
- **Idrofobici:** PLIP conta carboni di residui polari (Gln, Asp, Glu, Arg) e dell'eme; ligint si limita ai residui idrofobici.
- **Ponti d'acqua:** PLIP applica filtri propri (un donatore già in un H-bond diretto non forma un ponte; al massimo due ponti da
  donatore per acqua); ligint no, e riporta molti più ponti (precisione ≈ 40%).
- **Protonazione del ligando:** ligint assume le ammine alifatiche protonate e i carbossilati anionici; PLIP usa Open Babel.

## Ripetere la misura

```
pip install plip openbabel-wheel
python3 tests/benchmark.py cartella_con_pdb/ --groups gruppi.txt --out risultati/ [--all-copies] ["--extra=--plip-like"]
```

Scrive `report.md`, `disagreements.csv` (ogni disaccordo con atomi, distanze e causa probabile automatica) e
`per_category.csv` (TP/FP/FN per sito e tipo, per bootstrap e confronti appaiati). I report di questo repository sono in
`docs/benchmark/`: `AB`, `C`, `D_default`, `D_hb-long_salt-centroid`. Le cause probabili sono una classificazione
automatica di comodo (distanza oltre soglia, partner istidina...).

## Limiti

- Un solo strumento di riferimento; nessuna verità sperimentale.
- Nessun confronto per DNA/RNA (PLIP li ignora), alogeni con pochi casi (ordine di 10), cation-π pochi casi.
- Gli ID sono stati scelti a memoria da una sola persona: la rappresentatività dei gruppi non è stata verificata.
