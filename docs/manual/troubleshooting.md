# Troubleshooting

## The file does not open

**`Cannot open FILE`.** The path is wrong or the file cannot be read. Put the path in quotes if it contains spaces.

**`FILE is compressed (gzip): decompress it first, e.g. with gunzip`.** Files from the PDB archive are often `.pdb.gz` or `.cif.gz`.
Decompress them (`gunzip 1abc.cif.gz`) and use the result.

**The listing is empty, or a warning says that no atoms were read.** The file is not a structure ligint understands, or its format was
recognized wrongly. A PDBx/mmCIF file must contain an `_atom_site` category. You can force the format with `--format pdb` or `--format cif`.
See [Input formats](formats.md).

## `No ligand 'X' found`

- The name is the **residue name** of a `HETATM` record, in capital letters (`ATP`, not `atp`). Run `ligint file` without a name to list
  the ligands that the file really has.
- A **peptide** or a nucleic acid that binds the protein is part of the polymer (`ATOM` records) and is not a ligand for ligint.
- With a chain or a residue number the copy must exist: `ligint 1abc.pdb ATP B 501`. A residue number followed by `A` selects an
  insertion code, not a chain.
- The program exits with code 0 in this case: check the output in scripts.

## No interactions, or fewer than expected

- Look at the **Mean B-factor** line: a ligand with a B-factor much higher than the rest of the structure is poorly defined and its
  contacts are unreliable.
- If you chose a copy of the ligand that sits at the edge of the crystal, it may have few neighbors. Try another copy.
- Water bridges need `--no-water` to be absent, and appear only when the structure contains waters.
- Contacts beyond the thresholds are not listed: see the [thresholds](cli.md#thresholds). `--plip-like` or `--hb-long 4.1` accept longer
  hydrogen bonds when the geometry allows them.
- With several conformations of a residue, only one is used: check **Alternate conformation** (GUI) or `--altloc`.

## An expected hydrogen bond or salt bridge is missing

The roles of the ligand atoms are inferred from the geometry, and some rules are deliberate:

- An amide nitrogen is **not an acceptor**; a nitrogen with a lone pair in an aromatic ring is.
- An O or N bound to a **metal** is considered deprotonated and does not donate a hydrogen bond (`--metal-donors` changes it). A histidine
  bound to a zinc cannot make hydrogen bonds with the ligand.
- **Histidine is never charged**, so it makes no salt bridge.
- Aliphatic amines are assumed **protonated** and carboxylates **deprotonated**; amidines and guanidines not in a ring are cationic.

If the assumption is wrong for your ligand, give the exact chemistry with an SDF or MOL2 file (`--sdf`, `--mol2`, or **Add template...**).
If the template is **not applied**, the **Log** and the standard error say so: the file must have the same atoms at the same coordinates
as the structure, so use a file prepared from that complex.

## Many more water bridges than another program reports

ligint does not apply the extra filters that PLIP applies to water bridges (a donor already in a direct hydrogen bond cannot form a
bridge, at most two bridges per donor water). See [Comparison with PLIP](benchmark.md).

## Halogen bonds with fluorine are missing

Fluorine is not considered a halogen-bond donor in ligint (PLIP counts it). Chlorine, bromine and iodine are.

## A warning about invalid coordinates

`Warning: N atoms with invalid coordinates ignored`: some atoms have missing or absurd coordinates (for example `?` in an mmCIF file or
a value beyond 100,000 A) and were skipped.

## RCSB downloads and KLIFS pockets

**`RCSB: RCSB requires a four-character alphanumeric PDB ID`.** The ID is four letters or digits (`3EYG`); extended identifiers are not supported by the download.

**`RCSB HTTPS unavailable: rebuild with libcurl`.** This build of ligint has no network support: it can only open structures that are already in the cache or
that you downloaded yourself. Packages and builds with libcurl have it ([Downloading from the RCSB](rcsb.md#build)).

**`RCSB download failed ... (HTTP 0, Couldn't resolve host name)`.** No network connection or a proxy problem: libcurl uses the usual proxy variables
(`https_proxy`). A file that was already in the cache is not touched.

**`Use --download-pdb without a separate structure path`.** With `--download-pdb` give only the ligand (`ligint --download-pdb 3EYG STI`); the structure is the downloaded file.

**`KLIFS: Duplicate KLIFS position`, `Invalid KLIFS pocket position`, `Invalid PDB residue`, `KLIFS region missing`.** The mapping file is not valid: see
[the mapping file](klifs.md#the-mapping-file). No pocket file is written when there is an error.

**`KLIFS: KLIFS request failed` or `Online KLIFS unavailable`.** The service could not be reached, or this build has no network support. Use a mapping file
(`--klifs-mapping`); the online path of the native client is not verified ([KLIFS kinase pockets](klifs.md#the-klifs-service-not-verified)).

**The pocket is built but the figure has fewer than 85 cells.** The residues that are missing in the crystal (`_` in the mapping) are left out of the pocket. Check also that the chain and the residue numbers
of the mapping are the ones of the structure you analyze: a position whose residue is not in the structure is drawn as an empty cell.

## The graphical interface

**It does not start on Linux, with a message about the Qt platform plugin.** Install the Qt platform plugins (`qt6-qpa-plugins` on Debian and
Ubuntu); if the message mentions `xcb-cursor`, install `libxcb-cursor0`. To start the program only for tests, without a screen, use
`QT_QPA_PLATFORM=offscreen`.

**The fingerprint tab shows only a table.** This build of the program was made without the Qt Svg module. Everything else works, and the
SVG export still produces the figure.

**Windows shows "Windows protected your PC".** The installer is not code-signed. Choose **More info**, then **Run anyway**.

**Where are the settings?** The program remembers only the last folder you used, in the Qt settings of your system
(`~/.config/ligint/ligint-gui.conf` on Linux, the registry key `HKEY_CURRENT_USER\Software\ligint` on Windows).

## Reporting a problem

Open an issue in the [project repository](https://github.com/maxbar64/LigInt/issues) with: the output of `ligint --version`, the exact command (in
the GUI, copy it with **File > Copy command line**), what you expected and what happened, and, if you can share it, the structure file (or its PDB code).
For the graphical interface add the contents of the **Log** tab.
