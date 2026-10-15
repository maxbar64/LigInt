# Test data

Structures from the Protein Data Bank (public domain, CC0), used by the tests.

- `*.pdb`: files in the legacy PDB format, as downloaded from the PDB.
- `1AZM.cif`: the same entry as `1AZM.pdb`, **converted to PDBx/mmCIF with [gemmi](https://gemmi.readthedocs.io)**
  (`read_structure` + `make_mmcif_document`). It is not the file distributed by the PDB, but it has the usual structure of an
  mmCIF coordinate file (quoted values, text fields, `label_*` and `auth_*` columns, entities). It is used to check that the
  mmCIF reader gives the same interactions as the PDB reader.
- `1efy.sdf`: ligand of 1EFY with explicit hydrogens, generated with Protoss (ProteinsPlus).
