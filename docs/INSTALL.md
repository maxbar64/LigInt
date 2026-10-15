# Installation

Ready-made installers are built by the CI workflow (`.github/workflows/ci.yml`): they are attached to each
[GitHub release](https://github.com/maxbar64/LigInt/releases) and kept as artifacts of every workflow run.

| System | File | What it is |
|---|---|---|
| Linux (Debian, Ubuntu) | `ligint_<version>_amd64.deb` | programs, manual pages, menu entry and icon |
| Linux (Debian, Ubuntu) | `ligint-dev_<version>_amd64.deb` | headers, static library and CMake package, for programs that use the library |
| Linux | `ligint-<version>-linux-x86_64.tar.gz` | portable archive, no installation needed |
| Rocky Linux, AlmaLinux, RHEL 9 | `ligint-<version>-1.x86_64.rpm` | programs, manual pages, menu entry and icon |
| Rocky Linux, AlmaLinux, RHEL 9 | `ligint-devel-<version>-1.x86_64.rpm` | headers, static library and CMake package |
| Windows 10/11 (64 bit) | `ligint-<version>-win64.exe` | installer (wizard, Start menu, uninstaller) |
| Windows 10/11 (64 bit) | `ligint-<version>-win64.zip` | portable archive, no installation needed |

## Linux

```
sudo apt install ./ligint_0.8.2_amd64.deb          # the dependencies (Qt 6 libraries) are installed automatically
ligint --version
ligint-gui                                         # or from the applications menu: "ligint"
man ligint
sudo apt remove ligint                             # uninstall
```

- The `.deb` files are built on Ubuntu 24.04 (Qt 6.4) and depend on the Qt 6 packages of that release. On a distribution
  with different Qt packages build from source (see the README) or use `cpack` as shown below.
- `ligint-gui` is associated with PDB files (MIME type `chemical/x-pdb`): *Open with* shows it.
- The portable archive is extracted anywhere (`bin/ligint`, `bin/ligint-gui`); the Qt 6 libraries must be installed
  (`libqt6widgets6`, `libqt6svg6`, `qt6-qpa-plugins` on Debian/Ubuntu).

## Rocky Linux, AlmaLinux, RHEL 9

```
sudo dnf install epel-release                     # Qt 6 comes from EPEL on version 9
sudo dnf install ./ligint-0.8.2-1.x86_64.rpm      # the dependencies are installed automatically
ligint --version
sudo dnf remove ligint                            # uninstall
```

- The RPM is built by the CI in a Rocky Linux 9 container and installed on a clean AlmaLinux 9 container, so it follows the
  glibc and Qt of that release. It is **not** meant for version 8 or 10, and the `.deb` and the portable `.tar.gz` (built on
  Ubuntu 24.04) will not run on version 9 (they need a newer glibc).
- If Qt 6 is not available for the release, the CI builds the command line only: the RPM then has no `ligint-gui`
  (the workflow run shows a warning and says so in its summary).
- RHEL needs an active subscription and the CodeReady Builder repository (`crb`) to build from source.
- Other RPM-based systems (Fedora, openSUSE): build the RPM yourself with `cpack -G RPM` after building from source.

## Windows

1. Run `ligint-<version>-win64.exe` and follow the wizard: install folder, whether to add ligint to the `PATH`
   (so that `ligint` works in a terminal) and the Start menu entry.
2. Start **ligint** from the Start menu, or `ligint file.pdb` in a terminal.
3. To uninstall use *Settings > Apps*, or the *Uninstall* entry in the Start menu.

- **Silent installation:** `ligint-<version>-win64.exe /S` (optionally `/D=C:\Tools\ligint`, last on the command line).
- The installer is **not code-signed**: Windows SmartScreen may show "Windows protected your PC" the first time
  (*More info* > *Run anyway*). Signing needs a certificate, which this project does not have.
- The portable `.zip` can be extracted anywhere; run `bin\ligint-gui.exe` or `bin\ligint.exe`.
- The installer includes the Qt libraries and the Visual C++ runtime that the GUI needs.

## Building the installers yourself

The installers are produced by CPack from the same sources:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd build && cpack                 # Linux: .deb (+ ligint-dev) and .tar.gz;  Windows: NSIS installer and .zip
```

- Windows: use a "Developer Command Prompt" (MSVC), point CMake to Qt 6 with `-DCMAKE_PREFIX_PATH=C:\Qt\6.7.0\msvc2019_64`
  and install [NSIS](https://nsis.sourceforge.io) for the installer (the `.zip` does not need it). Qt 6.3 or later is required to bundle the Qt libraries.
- One generator only: `cpack -G DEB`, `cpack -G TGZ`, `cpack -G RPM`, `cpack -G NSIS`, `cpack -G ZIP`.
- Maintainer of the `.deb`/`.rpm`: `-DLIGINT_PACKAGE_MAINTAINER="Your Name <you@example.org>"` (the default is a placeholder);
  home page in the metadata: `-DLIGINT_HOMEPAGE=...`.
- **Network features** (the [RCSB download](RCSB.md) and the [KLIFS service](KLIFS.md)) are built when the development packages of libcurl (RCSB and KLIFS) and
  json-c (KLIFS) are found: `libcurl4-openssl-dev libjson-c-dev` on Debian and Ubuntu, `libcurl-devel json-c-devel` on RHEL and derivatives. The Linux packages built by the CI
  include them; the Windows installer does not (see [KLIFS kinase pockets](KLIFS.md) for a build with vcpkg). A build without them works, and still opens the
  structures that are in the download cache. The `ligint-dev` package then depends on those development packages, because its static library needs them.
- Without Qt only the command line is built and packaged (`-DLIGINT_BUILD_GUI=OFF`): the `.deb` then depends only on the C and C++ libraries.

## Publishing a release

Push a tag that starts with `v`; the workflow builds, tries and publishes the installers:

```
git tag v0.8.2 && git push origin v0.8.2
```

The version in the file names comes from `project(ligint VERSION ...)` in `CMakeLists.txt`: change it (and `LIGINT_VERSION`
in `include/ligint/ligint.hpp`) before tagging.

## What has been verified

| | Status |
|---|---|
| `.deb` (ligint, ligint-dev) | installed, run (CLI and GUI), a program built against `ligint-dev`, removed; `lintian` reports no errors or warnings |
| `.tar.gz` (Linux) | extracted in another folder and run |
| `.rpm` (Rocky/Alma/RHEL 9) | the layout (`lib64`, no system directories) was inspected here; **built in a Rocky 9 container and installed on AlmaLinux 9 by the CI, which has not run yet**; RHEL itself and versions 8 and 10 not tried |
| Windows installer (NSIS) | the installer script is generated and compiled (checked on Linux with a Linux payload); **not run on Windows**: CI builds and tries it there (silent install, run, uninstall) |
| Windows `.zip`, Qt bundling (`windeployqt`), MSVC build | **not tried by the authors**: only the CI run will tell |
| macOS | no installer (the library and the command line build; see the CI) |
