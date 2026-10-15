# KLIFS kinase pockets (experimental)

The [fingerprint](FINGERPRINT.md) of a kinase ligand is drawn on a **pocket definition**: a text file that says which residues form the pocket,
in which order and in which regions. [KLIFS](https://klifs.net) defines such a pocket for all the kinases, with 85 aligned positions. ligint can write its
pocket-definition file from a **KLIFS mapping** (the table of the 85 positions with the residues of a structure), without Python.

The file is a **structural definition** of the aligned binding pocket, not an annotation of ligand-protein contacts: the interactions are always computed by
ligint from the structure.

## The workflow

1. Get the KLIFS **mapping** of your structure (one row per position, see below): from a KLIFS export that you have, or from the service.
2. Generate the pocket-definition file:

   ```text
   ligint --klifs-pdb 3EYG --klifs-chain A --klifs-mapping KLIFS_pocket_export.tsv --klifs-output 3EYG-klifs-pocket.tsv
   ```

3. Draw the fingerprint of the ligand on it, with the **same structure and the same alternate conformation** as KLIFS:

   ```text
   ligint 3EYG.pdb LIG A 1 --altloc B --fingerprint svg --pocket 3EYG-klifs-pocket.tsv > 3EYG.svg
   ```

   KLIFS describes one *alternate model* of a chain (A, B, ...); in the PDB file that is normally an alternate location (altloc) of the atoms: choose it with
   `--altloc` (the listing of `ligint 3EYG.pdb` shows which ones the file has). Both steps can be one command: `ligint 3EYG.pdb LIG A 1 --altloc B
   --fingerprint svg --klifs-pdb 3EYG --klifs-chain A --klifs-mapping KLIFS_pocket_export.tsv`.

## The mapping file

A CSV or TSV file (the separator is recognized; quoted values, a byte order mark and Windows line ends are accepted; lines starting with `#` are
ignored) with the columns:

| Column | Meaning |
|---|---|
| `residue.klifs_id` | the KLIFS position, 1 to 85 (also accepted: `klifs_id`, `klifs_position`, `position`) |
| `residue.id` | the residue number of the structure, with its insertion code if any (`221A`); also `residue_id`, `pdb_residue`, `resi` |
| `residue.klifs_region` | the KLIFS region (`hinge`, `aC`, ...); also `klifs_region`, `region` |
| `residue.aa` | optional: the residue, in one letter or three (also `residue.amino_acid`, `residue_name`, `resname`) |

- `residue.id` must be the **author numbering** of the structure and chain that you analyze (the numbers of the PDB file), never an index along the
  sequence. A mmCIF file is read with the author numbers too ([Input formats](formats.md)).
- A residue marked `_` (or empty, `null`, `None`) is **missing in the crystal**: the position is left out of the pocket, never filled with an invented number.
  The figure then has fewer than 85 cells.
- A position outside 1 to 85, a position that appears twice, a residue number that is not an integer with an optional letter, or a missing region is an
  **error** and no file is written.
- Columns may come in any order and extra columns are ignored.

### A KLIFS `pocket.mol2`

`--klifs-mapping` also accepts the `pocket.mol2` file that KLIFS provides for each structure (it is recognized by its content). Its residues are, in order, the
positions 1 to 85 **that the structure has**; the regions are those of the KLIFS alignment, which ligint knows. A pocket with all 85 residues needs nothing more.
When some are missing in the crystal the file has fewer than 85 and does not say which, so the **pocket sequence** is needed: 85 characters, the one-letter code of each
position and `_` where the residue is missing (the field `pocket` of the KLIFS service and of the KLIFS overview table). Give it with `--klifs-sequence`. ligint
checks it against the file: the number of residues must agree, and every standard residue must be the one that the sequence has in that position, so a sequence that is
shifted or belongs to another structure is refused instead of producing a wrong pocket. The residue numbers are those of the file (the author numbering of the structure).

The generated file has the [pocket definition format](FINGERPRINT.md#pocket-definition-file): `@name`, one `@color` for each region (the colors of the
KLIFS regions) and one line per position with the chain, the residue, the position as label, the region and the residue name.

## Command line

| Option | Meaning |
|---|---|
| `--klifs-pdb ID` | the PDB ID (four letters or digits); only used for the name of the pocket when a mapping file is given |
| `--klifs-chain C` | the chain of the structure (**required**) |
| `--klifs-mapping FILE` | the mapping file or the `pocket.mol2` of KLIFS; without it the service is asked (see below) |
| `--klifs-sequence S` | the 85-character pocket sequence (`_` for a missing residue), for a `pocket.mol2` with fewer than 85 residues |
| `--klifs-output FILE` | where to write the pocket-definition file |
| `--klifs-model M` | the KLIFS alternate model, to choose among several (service only) |
| `--klifs-id N` | the KLIFS structure ID, to choose exactly (service only) |

Rules: with a structure file the generated pocket is used by `--fingerprint` (without it a warning says that the pocket is not used); with no structure
file `--klifs-output` is required; `--pocket` and `--klifs-*` cannot be given together. The messages begin with `KLIFS:` and the exit code of an error is 1.

## Graphical interface

In the **Pocket definition (fingerprint)** group, **Generate from KLIFS...** asks for the PDB ID (suggested from the file name), the chain and the source of the
mapping: **Local KLIFS mapping CSV/TSV** (a file chooser follows) or, **only in a build that has network support**, **Online (KLIFS service)** (then the KLIFS
model). It asks where to save the file and works in the **background**; when it finishes the pocket is loaded and the **Log** says how many positions were
mapped. A structure must be open first.

## The KLIFS service (partly verified)

**Verified by two real runs.** The list of the structures of a PDB entry works (two structures of one chain, IDs 1445 and 1453, were found and read), and the request for
the pocket answers with a **MOL2 file** (`@<TRIPOS>MOLECULE`, the name `HUMAN/JAK1_3eyg_altB_chainA`, 85 substructures), not with a table: ligint reads it as described above.
The reading was checked on four real KLIFS pocket files (BMX, BRAF and the two chains of ABL1) with the sequences of the KLIFS overview: the number of missing residues
agrees with what KLIFS says in all of them, and in ABL1 the hinge falls on Glu316-Phe317-Met318 with Thr315 as gatekeeper, as it should.

**Not yet confirmed**: that the answer about the structures names the pocket sequence `pocket` (ligint also tries `structure.pocket`) and the alternate model `alt` (and the other
names listed below), and a complete run to the end on a structure that has missing residues. If the sequence is not found the message lists the fields that the answer has;
with `--klifs-sequence` the command line can go on anyway.

Without a mapping file ligint asks the KLIFS web service, first for the structures of the PDB entry and then for the pocket of the one that matches the chain.

**What a real run showed.** The first request works against the live service: a user run on a real entry got the list of its structures and read the structure
ID and the chain from it (two structures of chain A, IDs 1445 and 1453). **Not confirmed**: the field that names the *alternate model* (ligint looks for `alt`,
`alternate_model`, `structure.alternate_model`, `alt_model` and `model`), and the second request, the one that returns the pocket (`structure_get_pocket`), including
the form of its answer. The address that the code uses (`https://klifs.net/api_v2/...`) is not one that I could find in the public documentation, where the examples
use `https://klifs.net/api/`; it answered, but treat the rest as unverified.

**Several structures for the same chain.** A chain can have one KLIFS structure for each alternate model. If the request does not say which, ligint stops and lists the
candidates with their IDs (`KLIFS matching structures: 2 (must be exactly 1). Candidates: [ID=1445 chain=A model=B] ...`):
choose with `--klifs-model M`, or with `--klifs-id N` when the answer does not name the model. In the graphical interface a window lists the candidates and lets you
pick one. KLIFS shows the alternate model of each ID on the page of the structure.

**When something does not match what ligint expects**, the message tells what the answer contained: the names of its fields, or the first characters of a
response that is not JSON. Please report them: they are what is needed to fix the reading.

If the request or the reading of the answer fails, the program says so and writes nothing; use a mapping file. The [Python helper](https://github.com/maxbar64/LigInt/blob/main/tools/README-KLIFS.md)
`tools/klifs_to_pocket.py` uses the OpenCADD client instead and is the better choice for the online step until the native one is verified.

Requests use HTTPS only, at most 3 redirects, a 10 s connection and a 25 s total timeout, and the answer is limited to 8 MB.

## Verified

By the automated tests, which use local files only: the reading of a `pocket.mol2` (all 85 residues, with a sequence, insertion codes, negative numbers, modified residues,
Windows line ends, every kind of error, random corruption, and the regions and colors, with their Greek letters, through the fingerprint and the SVG figure); the positions and the regions; missing residues; insertion codes; the residue names from the one-letter
code (all twenty); the same pocket from TSV, CSV, Windows line ends, byte order mark and permuted columns; the errors above (the file is not written);
corrupted files; the round trip through the pocket-definition format and the fingerprint (the positions, the regions and a hydrogen bond in the right cell);
the rules of the command line; the background task of the interface.

## Building with network support on Windows (vcpkg)

The service needs libcurl and json-c. From an x64 Native Tools command prompt, with vcpkg installed and its directory in `VCPKG_ROOT`:

```powershell
& "$env:VCPKG_ROOT\vcpkg.exe" install curl:x64-windows json-c:x64-windows
cmake -S . -B build-klifs -G Ninja `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DCMAKE_PREFIX_PATH="C:/Qt/6.9.3/msvc2022_64" `
  -DLIGINT_KLIFS_REQUIRE_ONLINE=ON
cmake --build build-klifs --config Release
```

Use a new build directory so that CMake finds the vcpkg toolchain at the first configuration. `LIGINT_KLIFS_REQUIRE_ONLINE=ON` makes a missing libcurl or
json-c an error instead of building without the service. With a dynamic triplet the DLLs must be available at run time. A successful build says nothing about
the compatibility with the live service.
