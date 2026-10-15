# Native KLIFS kinase pocket generator (experimental)

LigInt can generate its normal pocket-definition TSV from a KLIFS mapping without Python.
The TSV is a **structural definition** of the KLIFS aligned binding pocket (1–85),
not a ligand–protein contact annotation. Binding interactions are still computed by LigInt.

## Build

The offline importer uses only C++17. The online provider is enabled automatically when
`libcurl` and `json-c` development packages are discoverable using pkg-config.
If they are not installed, all other features continue to build and the KLIFS CSV/TSV
importer remains available.

## CLI (fully testable offline)

```
ligint --klifs-pdb 3EYG --klifs-chain A \
  --klifs-mapping KLIFS_pocket_export.tsv \
  --klifs-output 3EYG-klifs-pocket.tsv
ligint 3EYG.pdb --pocket 3EYG-klifs-pocket.tsv --fingerprint svg > fingerprint.svg
```

The input mapping needs columns `residue.klifs_id`, `residue.id`, and
`residue.klifs_region` (optionally `residue.aa`). Missing residues marked `_`
are omitted; numeric positions must be within 1..85 and duplicate mapped positions
are rejected. `residue.id` must use **PDB author numbering**, matching the input
protein chain. Never substitute contiguous sequence indices for author residue IDs.

## Online (experimental; live service compatibility must be confirmed)

```
ligint --klifs-pdb 3EYG --klifs-chain A --klifs-model B \
  --klifs-output 3EYG-klifs-pocket.tsv
```

The online implementation queries KLIFS `/api_v2` for a PDB entry and an associated
structure ID, then attempts to parse the pocket mapping response. **This path is not
verified against live KLIFS in the build environment**, where remote DNS is blocked.
KLIFS APIs may deliver PDB coordinate content for `structure_get_pocket` rather than
residue mapping records. In that case the generator fails explicitly rather than
manufacturing a mapping. Until the live response schema is verified, prefer the
local mapping export.

For ambiguous PDB/chain matches provide `--klifs-model` or `--klifs-id N`.
The CLI also accepts `--klifs-mapping`, and the GUI's `Generate from KLIFS…`
button supports both online and local CSV/TSV mapping sources. GUI online requests
currently block while waiting for the timeout; future releases should move them to
a worker.

## Validation

The C++ tests cover 1–85 position assignment, omitted residues, author numbering,
file roundtrip and duplicate KLIFS position rejection. No automated live network test
is included. GUI integration must be built and smoke-tested in a Qt-enabled environment.

## Windows: C++ online support using vcpkg (no Python)

From an x64 Native Tools command prompt, with vcpkg installed and its directory
recorded in `VCPKG_ROOT`:

```powershell
& "$env:VCPKG_ROOT\vcpkg.exe" install curl:x64-windows json-c:x64-windows
cmake -S . -B build-klifs -G Ninja `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DCMAKE_PREFIX_PATH="C:/Qt/6.9.3/msvc2022_64" `
  -DLIGINT_KLIFS_REQUIRE_ONLINE=ON
cmake --build build-klifs --config Release
```

Use a new build directory so CMake finds the vcpkg toolchain at first configure.
The option `LIGINT_KLIFS_REQUIRE_ONLINE=ON` intentionally fails configuration
if libcurl or json-c is missing rather than building an offline-only GUI.
For a dynamic vcpkg triplet, DLLs must be available at application runtime.

This resolves only dependency detection. The live KLIFS API response mapping
requires a separate end-to-end test; no guarantees about its compatibility are
made by a successful build.
