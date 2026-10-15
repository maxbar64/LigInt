# Input formats

ligint reads the structure from a **PDB** file or from a **PDBx/mmCIF** file. The format is recognized from the content
(an mmCIF file starts with `data_` or with a comment; a PDB file with a record name), not from the extension; `--format pdb`
or `--format cif` forces it. In the library: `Options::format` (`StructureFormat::Auto`, `Pdb`, `Mmcif`).

```
ligint 1abc.pdb LIG
ligint 1abc.cif LIG                 the same analysis, the same output
ligint structure.dat LIG --format cif
```

Why mmCIF: the legacy PDB format cannot represent the entries that the wwPDB is going to release (extended 12-character IDs,
five-character component ids for new ligands, more than 62 chains or 99,999 atoms). The wwPDB plans to stop generating PDB
files for new entries around July 2027, so mmCIF is the format to use for them. Results are the same for the same structure:
the analysis does not know which format was read.

## What is read from an mmCIF file

Only the atoms: the `_atom_site` category of the **first data block that has one**, and of the **first model**
(`pdbx_PDB_model_num`). The rest of the file is skipped without being interpreted.

| Item of `_atom_site` | Used for | If it is missing or null (`?`, `.`) |
|---|---|---|
| `Cartn_x`, `Cartn_y`, `Cartn_z` | coordinates | required (no atoms are read without them) |
| `auth_atom_id`, else `label_atom_id` | atom name | one of the two is required |
| `auth_comp_id`, else `label_comp_id` | residue name (up to 5 characters) | one of the two is required |
| `auth_asym_id`, else `label_asym_id` | chain (can have several characters) | one of the two is required |
| `auth_seq_id`, else `label_seq_id` | residue number | 0 |
| `pdbx_PDB_ins_code` | insertion code | none |
| `label_alt_id` | alternate conformation (`--altloc`) | none |
| `type_symbol` | element | first letter of the atom name |
| `group_PDB` | `HETATM` marks the ligands and the other non-polymer residues | all atoms are considered `ATOM` |
| `B_iso_or_equiv` | B-factor of the ligand (warning when it is much higher than the rest of the structure) | 0 |
| `pdbx_PDB_model_num` | first model only | single model |

The **author** identifiers (`auth_*`) are used first because they are the ones of the PDB file, of the literature and of
databases such as KLIFS; the `label_*` ones are the fallback. Tags are case-insensitive, in any order. The file is read as a
stream: structures with millions of atoms are not held in memory as text.

Not read: bonds and bond orders (`_chem_comp_bond` is a dictionary, not part of the coordinate file), formal charges
(`pdbx_formal_charge`), occupancies, anisotropic factors, `_struct_conn`. As with PDB files, the chemistry of the ligand is
inferred from the geometry unless an SDF or MOL2 template is given.

## Residue labels with long chain names

A chain of one character keeps the usual label, `ASP 189A` (name, number, insertion code if any, chain). A chain of several
characters, which exists only in mmCIF, is written in brackets so that the label stays unambiguous: `ASP 189 [AA]`,
`ASP 189B [AA]` (with insertion code B). The same rule is used in the interaction tables, in the CSV, in the fingerprints and in
the ligand ids (`LIG 1 [LL]`). The chain given on the command line (`ligint x.cif LIG LL 1`) and the chain of a pocket definition
file (`LL  1  p1  reg`) are the plain chain names.

## Limits

- **Compressed files are refused.** Files from the PDB archive are often `.cif.gz` or `.pdb.gz`: ligint stops with a message
  (`... is compressed (gzip): decompress it first, e.g. with gunzip`) instead of reading garbage. There is no zlib dependency.
- **mmCIF files without `_atom_site`** (for example a dictionary or a structure-factor file) give a warning and no atoms.
- **BinaryCIF**, mmTF and other formats are not supported.
- **What the reader was checked on.** (1) 119 real structures converted from PDB to mmCIF with gemmi: 1314 runs, results identical to the PDB
  ones. (2) Seven model files written by the wwPDB, taken from its
  [extended-identifier examples](https://github.com/wwPDB/extended-wwPDB-identifier-examples) (280 KB to 1.3 MB; five-character ligand ids, glycans,
  alternate conformations, a `data_pdb_...` block name): read, and the results are identical to those obtained from a PDB
  conversion made by gemmi (44 comparisons; the ligand names were shortened to three characters for the PDB side). (3) Synthetic and
  randomly corrupted files, and a 103 MB file with two million atoms (the reading takes about 3 s). It has **not** been tried on the large
  entries of the PDB archive (ribosomes, virus capsids) or on files with unusual `_atom_site` columns: if one does not read correctly, send it.
