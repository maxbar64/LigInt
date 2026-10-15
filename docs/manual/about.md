# About

**ligint** finds non-covalent ligand-protein interactions in PDB and PDBx/mmCIF files, as a command-line program, a graphical interface and
a C++ library. The version of this manual is shown at the bottom of every page and printed by `ligint --version`.

## License

ligint is distributed under the MIT license: see [LICENSE](../../LICENSE). The structures in the test data come from the Protein Data Bank,
whose data are in the public domain (CC0).

## Project

The source code, the issue tracker and the releases are at [github.com/maxbar64/LigInt](https://github.com/maxbar64/LigInt).
Contributions are welcome: see [CONTRIBUTING.md](../../CONTRIBUTING.md) and the page [For developers](developers.md).

## Acknowledgments and related tools

- [PLIP](https://github.com/pharmai/plip), the Protein-Ligand Interaction Profiler, was used as the reference for the comparison reported in
  [Comparison with PLIP](benchmark.md). ligint is an independent implementation, not a port of PLIP.
- [KLIFS](https://klifs.net) inspired the pocket-position fingerprints; ligint does not include KLIFS data.
- The graphical interface is built with [Qt](https://www.qt.io). [gemmi](https://gemmi.readthedocs.io) was used to generate mmCIF files for
  the tests.
- The Protein Data Bank ([rcsb.org](https://www.rcsb.org), [wwpdb.org](https://www.wwpdb.org)) provides the structures.

## The manual

The pages of this manual are written in Markdown (`docs/manual/` and the other files of `docs/`) and converted to HTML with
`tools/build_manual.py`. The HTML is embedded in the program (Help menu) and installed with it in the documentation folder
(`/usr/share/doc/ligint/manual/` on Linux), where you can open `index.html` with any web browser.
