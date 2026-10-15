# Native RCSB downloads (C++17 / libcurl)

The File menu and Structure panel contain **Download from RCSB…**. Enter a four-character PDB accession and select PDB or mmCIF. LigInt retains all downloaded ligands, ions and waters and loads the saved file into the normal structure browser. KLIFS annotations are independent and optional.

## CLI

```sh
ligint --download-pdb 3EYG --download-only
ligint --download-pdb 3EYG --download-format cif --download-cache my_cache
ligint 3EYG.pdb STI --fingerprint svg
```

The last example should select the ligand through the ordinary positional arguments **after loading**, not by appending a positional argument to `--download-pdb`: to process a selected ligand use the downloaded file path as the structure argument in a second invocation.

Options: `--download-pdb ID`, `--download-format pdb|cif`, `--download-output FILE`, `--download-cache DIR`, `--download-refresh`, `--download-only`. Default cache is the current directory for the CLI and the application-local data directory in the GUI. The cache contains actual PDB/mmCIF files, not just metadata. Cache files are validated before reuse. Downloads are written to a temporary file, checked, and renamed on success. HTTPS/TLS certificates are verified by libcurl.

## Build

With a vcpkg curl port installed, configure CMake with the vcpkg toolchain, and `-DLIGINT_RCSB_REQUIRE_ONLINE=ON` so missing HTTP support is treated as an error. RCSB download only requires libcurl; JSON-C is needed separately for online KLIFS.

This experimental GUI path performs the transfer synchronously (it can briefly block input). The GUI Qt build and live network transfer have not been validated in the delivery environment.
