# GPCR pockets from generic residue numbers (GPCRdb)

For the [fingerprint](FINGERPRINT.md) and the [interaction matrix](MATRIX.md) of a G protein-coupled receptor (GPCR) the positions of the pocket are **generic residue
numbers**: the same number names the equivalent residue in every receptor of a class, whatever the residue number of the protein. ligint reads the numbers that
[GPCRdb](https://gpcrdb.org) puts in a structure, takes the pocket that you ask for, and writes the [pocket definition](FINGERPRINT.md#pocket-definition-file).

## The numbers

- **GPCRdb scheme**, written with an `x`: `3x32` is the position 32 of helix 3. It is built from crystal structures and corrects the helix bulges and constrictions, so
  that equivalent residues have the same number in all the receptors. It exists for every class.
- **Ballesteros-Weinstein (BW)**, written with a point: `3.32`. It is sequence-based and exists for class A (class B, C and F have their own schemes: Wootten, Pin, Wang).
- In both, the helix comes first and the most conserved residue of the helix is `.50` (`3x50`). A **bulge** residue repeats the number of the one before, followed by a `1`: `2x551`.
- The **segments**: the helices `TM1`-`TM7` (numbers `1x..` to `7x..`), helix 8 (`8x..`) and the loops, which have two digits: `12` is ICL1 (between TM1 and TM2), `23` ECL1,
  `34` ICL2, `45` ECL2, `56` ICL3, `67` ECL3. The segment is the region of the pocket definition, so the fingerprint and the matrix draw it.

**The two schemes are not the same near bulges and constrictions.** In the receptor that I used to check this (the human beta-2 adrenergic receptor, 2RH1, numbered by GPCRdb) the
GPCRdb number `5x42` is Ala202 and `7x38` is Asn312, `7x42` is Tyr316; in the BW numbering that appears in the literature the same residues are Ser203 (5.42), Asn312 (7.39) and Tyr316 (7.43)
(check them on GPCRdb for your receptor). A list of residues taken from a paper in BW numbers is **not** the same list in GPCRdb numbers, and ligint does not convert one into the other: it
reads the BW numbers only where GPCRdb wrote them (see below) and refuses to guess.

## The pocket is what you choose

A kinase has one pocket of 85 positions ([KLIFS](klifs.md)). A GPCR has numbers for its whole 7TM domain, and the pocket depends on the class and on the site (orthosteric, allosteric...). So the pocket
of interest is **a list of generic numbers that you give**, for example the residues of the site that the ligands of your series share. The pocket definition has those positions, ordered along the
receptor (TM1, ICL1, TM2, ECL1, TM3...), and each one with its segment as region.

## Getting the numbers of a structure

**From GPCRdb.** Its tool *Generic residue numbering (PDB)* ([upload a PDB file](https://gpcrdb.org/structure/generic_numbering_index)) annotates the file with the numbers in the B-factor column and returns it. ligint reads
the annotated file as the tool writes it:

| Where | What |
|---|---|
| the **CA** atom of a residue | the GPCRdb number, for example `3.32` for `3x32`; for the loops `34.50` is `34x50` |
| the **N** atom | the Ballesteros-Weinstein number, if the file has it (class A) |
| a B-factor of `0.00` | the residue has no number |
| a **negative** B-factor | a bulge: the number of the preceding residue followed by `1` (`-2.55` is `2x551`) |
| a number that repeats on the next residue | also a bulge (files that were rewritten by another program lose the negative sign): `2x55` on two residues in a row is `2x55` and `2x551` |
| a B-factor of `1.00` | a filler, not a number |

A file that does not carry numbers is **refused** (a PDB file with its normal B-factors has a few residues that look like numbers, but not many, and they do not grow along the helices): the message tells
which tool to use. The residues with a number must be at least 20 and 30% of the chain, and the numbers must grow along the helices.

**From a table.** A TSV or CSV with a header line, one line per residue (the columns can be in any order; other columns are ignored):

| Column | Meaning |
|---|---|
| `resi` | the residue number in the structure (`113`, `316A`); also `resseq`, `residue`, `resnum`, `number` |
| `gpcrdb` | the GPCRdb number (`3x32`); also `generic`, `generic_number`, `display_generic_number` |
| `bw` | the Ballesteros-Weinstein number (`3.32`) |
| `chain`, `segment`, `aa` | optional: the chain, the segment (`TM3`...; otherwise it comes from the number), the amino acid (one or three letters) |

A cell may have both numbers as GPCRdb writes them (`3.32x32`) and `-` or an empty cell means none. The residue numbers must be those of the structure that you analyze (a receptor with a fusion partner
or a different numbering than the table needs the table of that structure). Use this when you have the numbers from a lookup table of GPCRdb or from your own work.

## Command line

```text
ligint --gpcrdb-pdb numbered.pdb --gpcrdb-positions 3x32,3x33,3x36,6x48,7x39 --gpcrdb-output pocket.tsv
ligint --gpcrdb-mapping numbers.tsv --gpcrdb-chain A --gpcrdb-positions positions.txt --gpcrdb-output pocket.tsv
ligint 2RH1.pdb CAU --fingerprint svg --gpcrdb-pdb numbered.pdb --gpcrdb-positions positions.txt > figure.svg
ligint numbered.pdb CAU --fingerprint csv --gpcrdb-positions 3x32,3x33,3x36        the structure of the analysis is the numbered file
```

| Option | Meaning |
|---|---|
| `--gpcrdb-positions LIST` | **required**: the pocket, generic numbers separated by commas, blanks or lines (`3x32,3x33`), or the name of a file with them (`#` starts a comment) |
| `--gpcrdb-pdb FILE` | the structure annotated by GPCRdb; by default the structure of the analysis |
| `--gpcrdb-mapping FILE` | a table of numbers instead of the annotated structure |
| `--gpcrdb-chain C` | the chain of the receptor; by default the only chain that has numbers |
| `--gpcrdb-name TEXT` | the name of the pocket definition |
| `--gpcrdb-output FILE` | where to write the pocket definition (required without a structure) |

With a structure and `--fingerprint` the pocket is used at once; with no structure the file is written. A position that this structure does not have (a residue that is missing or not numbered) is listed as a warning and left out, so the same list works for
a series of structures; if none is there it is an error. `--pocket` and `--klifs-*` cannot be given with `--gpcrdb-*`.

For a **matrix** of several receptors write the pocket of each structure with the same list of positions and name it in the `pocket` column of the [list of ligands](matrix.md#the-list-of-ligands): the columns are the generic numbers.

## Example: carazolol in the beta-2 adrenergic receptor

`docs/examples/b2ar-positions.txt` has 18 GPCRdb positions around the site of the beta blockers (an illustration, not a curated definition). With the structure 2RH1 and a copy of the receptor annotated by GPCRdb:

```text
ligint 2RH1.pdb CAU --fingerprint csv --gpcrdb-pdb numbered.pdb --gpcrdb-positions docs/examples/b2ar-positions.txt
```

gives, among the 18 positions, a hydrogen bond and a salt bridge with Asp113 (`3x32`), hydrophobic contacts with Val114 (`3x33`), Val117 (`3x36`), Phe289 (`6x51`) and Phe290 (`6x52`), hydrogen bonds with Ser203 (`5x43`),
Asn312 (`7x38`) and Tyr316 (`7x42`), and an aromatic edge-to-face contact with Phe290.

## Checked and not checked

- **Checked with a real structure.** A GPCR with numbers in the B-factor (the beta-2 adrenergic receptor), taken from a public repository of a tool that uses GPCRdb numbers: the helices, the loops (12, 23, 34, 45), helix 8, the zeros for the residues without a number, the repeated numbers of the bulges
  (2x55 and 5x46 appear twice, so 2x551 and 5x461) and the known residues (Asp113 is `3x32`, Trp286 `6x48`, the DRY of `3x50`). With that file and 2RH1 the fingerprint of carazolol has the contacts above. The original 2RH1, without numbers, is refused.
- **Checked only on files that I built from the description of the tool**: the BW numbers in the N atoms, the negative B-factors of the bulges, the table in its various forms, several chains. The real file that I had
  did **not** have the BW numbers (the N atoms were `1.00`) and had no negative values, because it had been rewritten by another program: the output of the GPCRdb tool itself, as it is, has not been tried.
- **Not done**: downloading the numbers from the GPCRdb web service (the files are made with its tool and you give them to ligint), the interface (the **Matrix** and the pocket group of the window do not have a *Generate from GPCRdb* yet),
  classes B, C and F beyond what the GPCRdb scheme gives (their sequence-based numbers are read only if the tool writes them as `N.NN` in the N atoms), annotated mmCIF files (the reader accepts them but it was not tried).
