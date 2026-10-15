# Downloading structures from the RCSB

ligint can download a structure from the [RCSB Protein Data Bank](https://www.rcsb.org) and analyze it in the same step, in the graphical interface
and on the command line. The download needs **libcurl** at build time; a build without it can still use the structures that are already in the
cache (see [Installation](INSTALL.md)). The file is the ordinary PDB or PDBx/mmCIF file of the entry, with all its ligands, ions and waters.

## Graphical interface

**File > Download from RCSB...** (or the button of the same name in the **Structure** group) asks for the four-character PDB ID and the format (PDB or
mmCIF). The transfer runs in a **background thread**: the window stays responsive and a busy cursor shows that it is working. When it ends the
structure is opened; choose a ligand and press **Analyze**. The file is kept in a folder of your user data (the **Log** tab shows its path), so opening the
same entry again does not download it again. If the transfer fails, the message is shown and written in the **Log**, and nothing is opened.

## Command line

```text
ligint --download-pdb 3EYG                       download (or find in the cache) and list the ligands
ligint --download-pdb 3EYG STI --csv             download and analyze the ligand STI, in one command
ligint --download-pdb 3EYG STI A 201 --altloc B  chain, residue number and the other options work as usual
ligint --download-pdb 3EYG --download-only       only download; prints the path of the file
ligint --download-pdb 3EYG --download-format cif --download-cache my_cache
ligint --download-pdb 3EYG --download-output 3eyg.pdb --download-refresh
```

With `--download-pdb` the arguments after the options are the **ligand** (`RESNAME [chain [resseq]]`), because the structure is the downloaded
file; giving also the path of an existing structure file is an error.

| Option | Meaning |
|---|---|
| `--download-pdb ID` | the four-character PDB ID (letters or digits; lower case is accepted) |
| `--download-format pdb\|cif` | the file format (default `pdb`) |
| `--download-output FILE` | where to save the file (default `ID.pdb` or `ID.cif` in the cache folder) |
| `--download-cache DIR` | the cache folder (default: the current folder) |
| `--download-refresh` | download again even if a valid file is already there |
| `--download-only` | stop after the download and print the path of the file |

The path of the file is written to the standard error (`RCSB structure: ...`); with `--download-only` it is also written to the standard output.

## How it works

- The **ID is checked** before anything else: exactly four letters or digits, so no path or address can be built from it.
- The file is requested over **HTTPS only**, with at most 3 redirects (HTTPS only too), a 10 s connection timeout and a 90 s total timeout; the
  certificates are verified by libcurl.
- The data go to a **temporary file**, which is checked (a PDB file must contain a header or atom records, an mmCIF file a `data_` block) and only then
  replaces the destination. If the transfer fails or the data are not a structure, the temporary file is removed and a file that was already in the
  cache is **left untouched**.
- A file already in the cache is used if it is valid (not empty, at least 100 bytes, with the expected content); a damaged entry is never trusted.
- Compressed files are not downloaded: the RCSB service returns plain text for these URLs.

## Build

The download is built when libcurl is found (CMake looks for it with pkg-config or `find_package(CURL)`); `-DLIGINT_RCSB_REQUIRE_ONLINE=ON` makes a
missing libcurl an error instead of a warning. Both libcurl 7.85 and later and older versions (for example 7.76, the one of RHEL 9) are supported. JSON-C
is not needed for the downloads (it is for the [KLIFS service](KLIFS.md)). In a build without libcurl every `--download-*` request that is not
already in the cache fails with `RCSB HTTPS unavailable`.

## Verified and not verified

Checked by the automated tests, without any network traffic: the ID validation, the cache, the choice of the output file, the arguments of the
command, the error paths (a transfer sent to a proxy that does not exist leaves the cache and the folder unchanged), the background task of the
interface. **Not checked: a real transfer from the RCSB servers**, which the test environment cannot reach; try `ligint --download-pdb 1AZM --download-only`
and tell if it fails.
