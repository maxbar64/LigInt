# Command-line reference

```text
ligint structure [RESNAME [chain [resseq[A|-]]]] [options]
```

`structure` is a PDB or PDBx/mmCIF file (`.pdb`, `.ent`, `.cif`, `.mmcif`; the format is recognized from the content). Compressed
files (`.gz`) are refused: decompress them first.

## Arguments

| Argument | Meaning |
|---|---|
| `structure` | the file to read. Alone, it **lists the ligands** (`HETATM` residues) with their number of copies, and the alternate conformations of the file if it has any. |
| `RESNAME` | the name of the ligand to analyze (up to 5 characters in mmCIF), for example `ATP`. Every copy of the ligand in the file is analyzed. |
| `chain` | restricts the analysis to the copy in this chain (a chain of several characters is allowed in mmCIF). |
| `resseq` | restricts it to this residue number. A trailing letter selects an insertion code (`501A`); a trailing `-` selects the residue without insertion code (`501-`); the bare number matches any code. |

If the ligand is not found, `ligint` prints the list of the ligands present and `No ligand 'X' found.`, and **still exits with code 0**:
check the output, not only the exit code, in scripts.

## Options

| Option | Meaning |
|---|---|
| `--csv [OUTPUT.csv]` | write interactions as CSV to the named file, or to standard output without a filename; diagnostic information goes to standard error |
| `--sdf FILE`, `--mol2 FILE` | chemical annotation of the ligand (bond orders, explicit hydrogens, formal charges); the ligand is matched to the molecule by coordinates ([details](interactions.md#sdf-and-mol2)) |
| `--no-water` | ignore the water molecules (no water bridges) |
| `--no-merge` | do not merge covalently linked `HETATM` residues (glycans, dipeptides) into one ligand |
| `--plip-like` | preset for the comparison with PLIP: `--hb-long 4.1 --salt-centroid 5.5 --hal-max 4.0`; explicit options take precedence |
| `--metal-donors` | let O and N coordinated to a metal remain possible hydrogen-bond donors |
| `--backbone-hyd` | also count apolar contacts with the protein backbone |
| `--altloc X` | use the alternate conformation `X` for the residues that have it; the others keep their first conformation. The program reports for how many residues `X` was used and warns if the file has no conformation `X` |
| `--format F` | `auto` (default: recognized from the content), `pdb` or `cif` |
| `--fingerprint FORMAT`, `--pocket FILE` | write the [interaction fingerprint](fingerprint.md) (`csv`, `json` or `svg`), optionally with a pocket definition; can be combined with named CSV and text reports |
| `--download-pdb ID` | download the structure from the RCSB instead of giving a file: the arguments after the options are then the ligand. `--download-format pdb\|cif`, `--download-output FILE`, `--download-cache DIR`, `--download-refresh` and `--download-only` (stop after the download and print the path; it requires `--download-pdb`) complete it ([Downloading from the RCSB](rcsb.md)) |
| `--gpcrdb-positions LIST`, `--gpcrdb-pdb FILE`, `--gpcrdb-mapping FILE`, `--gpcrdb-chain C`, `--gpcrdb-name TEXT`, `--gpcrdb-output FILE` | build the pocket definition of a GPCR from generic residue numbers (the positions that you list), for `--fingerprint` or to save ([GPCR pockets](gpcrdb.md)) |
| `--matrix LIST` | the [interaction matrix](matrix.md) of the ligands of a list file instead of a structure: one row for each ligand, one column for each pocket position and class. `--matrix-format csv\|json\|svg` (default `csv`), `--matrix-classes klifs\|all\|NAMES`, `--matrix-hide-empty`, `--matrix-title TEXT` and `--pocket FILE` (the common pocket definition) complete it |
| `--klifs-pdb ID`, `--klifs-chain C`, `--klifs-mapping FILE`, `--klifs-sequence S`, `--klifs-output FILE`, `--klifs-model M`, `--klifs-id N` | build a pocket definition from a KLIFS mapping, to use with `--fingerprint` or to save ([KLIFS kinase pockets](klifs.md)) |
| `--version` | print the version |
| `--help`, `-h` | print the usage and the default values of the thresholds |

An unknown option, an option without its value, a non-numeric threshold, a missing file or a malformed residue number ends the program with
an error message and exit code **1**.

## Thresholds

Every criterion has a threshold with a default value that can be changed with `--name value`. `ligint --help` prints the defaults of your version.
Thresholds above 6 A, which the grid search cannot reach, are limited to 6 A with a warning. The meaning of each criterion is explained in
[Interactions and assumptions](interactions.md).

| Option | Default | Meaning | In the GUI |
|---|---|---|---|
| `--hb-min` | 2.2 | smallest donor-acceptor distance of a hydrogen bond (A) | |
| `--hb-max` | 3.5 | largest donor-acceptor distance (A) | **H-bond, D...A max** |
| `--hb-ha` | 2.7 | largest H...acceptor distance (A); applies when the hydrogens are known | |
| `--hb-angle` | 120 | smallest D-H...A angle (degrees); applies when the hydrogens are known | |
| `--hb-ante` | 90 | without known hydrogens: smallest angle (degrees) between a neighbor of the donor or acceptor, that atom and the partner | |
| `--hb-long` | 0 (off) | also accept donor-acceptor distances up to this value (A) when the geometry lets a hydrogen point at the acceptor | **Long H-bond, up to** |
| `--hb-long-tol` | 40 | angular tolerance (degrees) of the geometry test of long hydrogen bonds | |
| `--hyd-min` | 3.3 | smallest carbon-carbon distance of a hydrophobic contact (A) | |
| `--hyd-max` | 4 | largest carbon-carbon distance of a hydrophobic contact (A) | **Hydrophobic max** |
| `--salt-max` | 4 | largest distance between oppositely charged atoms of a salt bridge (A) | **Salt bridge (atoms)** |
| `--salt-centroid` | 0 (off) | judge salt bridges between the centers of the charged groups, up to this distance (A) | **Salt bridge (centroids)** |
| `--rep-max` | 4.5 | largest distance of an electrostatic repulsion between like charges (A) | |
| `--hal-max` | 3.5 | largest distance of a halogen bond (A) | |
| `--hal-angle` | 140 | smallest C-X...acceptor angle of a halogen bond (degrees) | |
| `--metal-max` | 3.5 | upper limit of a metal-coordination distance (A); the limit of each element (Mg 2.5, Zn 2.6, Ca 2.8...) applies within it | |
| `--pi-max` | 5.5 | largest distance between the centers of two aromatic rings (A) | |
| `--pi-offset` | 3.5 | largest lateral offset of the rings in parallel stacking (A) | |
| `--pi-angle` | 30 | largest angle between the normals of the rings in parallel stacking (degrees) | **pi-stacking, max angle** |
| `--cp-max` | 6 | largest distance between a cation and the center of a ring for cation-pi (A) | |
| `--clash` | 0.5 | a steric clash is a distance shorter than the sum of the van der Waals radii minus this value (A) | |
| `--wb-min` | 75 | smallest ligand-O(water)-protein angle of a water bridge (degrees) | |
| `--wb-max` | 140 | largest ligand-O(water)-protein angle of a water bridge (degrees) | |

## Output

### Text (default)

For each ligand the program prints the annotation and B-factor information, a table, two summaries and, when the ligand touches
metals, their coordination sphere (see the example in the [Quick start](quickstart.md)). The columns of the table are the **type** of the
interaction, **Fav** (`+` favorable, `-` clash or repulsion), the **ligand atom**, the **protein atom** (as `RES numberCHAIN:ATOM`), the
distance **d(A)** and **Notes**. A chain of several characters (mmCIF) is written in brackets: `ASP 189 [AA]:OD1`.

### Multiple output reports

The three reports are independent. Use `--text OUTPUT.txt` for the complete readable analysis, `--csv OUTPUT.csv` for detailed interaction records, and `--fingerprint OUTPUT.csv|OUTPUT.json|OUTPUT.svg` for the interaction fingerprint. They can be requested in a **single analysis**:

```bash
ligint receptor.pdb --mol2 ligand.mol2 \
  --text interactions.txt \
  --csv interactions.csv \
  --fingerprint fingerprint.json
```

With no report option, the readable text report is printed to stdout as before. Bare `--csv` continues to print CSV to stdout; bare `--fingerprint csv|json|svg` also prints to stdout when selected by itself. When several reports are requested, at most one may target stdout; the remaining reports need distinct file paths. `--text` always requires a filename. The optional filename following `--csv` must end in `.csv`.

### CSV (`--csv`)

One header line and one row per interaction, with these columns:

| Column | Content |
|---|---|
| `ligand` | the ligand copy, for example `AZM 262A` |
| `type` | the interaction type (the list is in [Interactions and assumptions](interactions.md)) |
| `favorable` | `1` or `0` |
| `ligand_atom`, `protein_atom` | the atoms, as in the text output |
| `distance` | in angstroms |
| `notes` | the notes of the text output |

Rows are sorted by type, ligand atom and distance. The information lines and the warnings go to standard error, so redirecting
standard output gives a clean file.

### Fingerprints

`--fingerprint csv|json|svg` selects the fingerprint stream: see [Interaction fingerprints](fingerprint.md).

## Exit status

| Code | Meaning |
|---|---|
| 0 | the analysis ran (also when no ligand matched: see above), or `--help` / `--version` |
| 1 | error: file that cannot be opened or is compressed, unknown option, missing or non-numeric value, malformed residue number, unknown `--format` or `--fingerprint` value, `--pocket` without `--fingerprint`, an invalid PDB ID or a failed download (`RCSB: ...`), an invalid KLIFS mapping or a missing `--klifs-chain` (`KLIFS: ...`) |

## Examples

```text
ligint 1abc.pdb                                       list the ligands
ligint 1abc.pdb ATP A 501 --csv > atp.csv             one copy of ATP, as CSV
ligint 1abc.pdb LIG --plip-like --no-water            the PLIP-like preset, without water bridges
ligint 1abc.pdb LIG --hb-max 3.2 --hyd-max 4.2        stricter hydrogen bonds, looser hydrophobic contacts
ligint 1abc.pdb LIG --sdf lig.sdf                     exact chemistry from an SDF of the ligand
ligint 3eyg.cif LIG A 1 --altloc B --fingerprint svg --pocket pocket.tsv > figure.svg
ligint --download-pdb 3EYG STI --csv                  download from the RCSB and analyze the ligand in one command
ligint --klifs-pdb 3EYG --klifs-chain A --klifs-mapping map.tsv --klifs-output 3eyg-pocket.tsv     a KLIFS pocket definition from a mapping
ligint --matrix ligands.tsv --pocket pocket.tsv --matrix-format svg > matrix.svg     the interaction matrix of the ligands of a list
ligint numbered.pdb CAU --fingerprint svg --gpcrdb-positions 3x32,3x33,6x48 > figure.svg     the fingerprint on a GPCR pocket, by generic residue numbers
```

Analyze many structures that contain the same ligand name, in a POSIX shell:

```text
for f in complexes/*.pdb; do ligint "$f" LIG --csv > "${f%.pdb}.csv"; done
```

Keep only the hydrogen bonds of a CSV: `ligint 1abc.pdb LIG --csv | awk -F, 'NR == 1 || $2 == "Hydrogen bond"'`.

To repeat an analysis made in the [graphical interface](gui.md#reproducing-an-analysis-on-the-command-line), copy the command line that it
writes in its log.

## External ligand poses (docking)

Batch summaries (`--batch-poses csv|json`) analyze every imported molecule independently and report per-pose favorable contact counts, hydrogen bonds, hydrophobic and aromatic contacts, salt and water bridges, clashes, repulsions, possible covalent contacts and the number of interacting residues. These are descriptive counts, **not a calibrated affinity or docking score**. A failed pose is recorded with an error and the process returns a nonzero status. Batch mode cannot be combined with `--fingerprint`, `--csv` or a specified crystallographic ligand.

When no ligand residue is supplied, `--sdf` or `--mol2` analyzes the **first** molecule by default (or the 1-based molecule selected using `--external-index N`) in the external file directly as a ligand, even if the receptor structure contains no matching HETATM residue. Coordinates, bonds, hydrogens and charges come from the external molecule; the receptor comes from the PDB/mmCIF file. The two inputs must use the same coordinate frame. When a ligand residue name is given, the existing coordinate-matched template behavior is unchanged.

```sh
ligint receptor.pdb --sdf docked_pose.sdf --csv
ligint receptor.pdb --mol2 docked_pose.mol2 --fingerprint json
ligint receptor.pdb --sdf multiple_poses.sdf --external-index 2 --fingerprint json
ligint receptor.pdb --sdf multiple_poses.sdf --batch-poses csv > pose_summary.csv
ligint receptor.pdb --sdf multiple_poses.sdf --batch-poses json > pose_summary.json
ligint complex_with_original_ligand.pdb --sdf docked_pose.sdf --exclude-pdb-ligand BZC --csv
```
