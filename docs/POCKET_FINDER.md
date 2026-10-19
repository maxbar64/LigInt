# Experimental buriedness pocket finder

`ligint-pocket` is a standalone CLI for ligand-independent geometric pocket discovery.
It is built with `LIGINT_BUILD_CLI=ON`. It does **not** modify the regular `ligint` interaction analysis.

```bash
ligint-pocket protein.pdb --chains A,B --output pockets.json --points voxels.csv --pdb pockets.pdb
ligint-pocket protein.cif --format cif --spacing 1.0 --directions 42 --top 10
```

The input is parsed using the project's PDB/mmCIF readers (first structural model).
The first version selects heavy atoms from `ATOM` records (not HETATM, including waters and
nonstandard amino acids). Chain IDs use the author's nomenclature; multi-character mmCIF
chain IDs are accepted.

## Algorithm (experimental)

1. Construct a union of atomic van der Waals volumes on a regular Cartesian grid.
2. Sample unoccupied grid points within a configurable distance of the protein.
3. Cast uniformly distributed rays and compute the fraction that intersect protein volume
   within the ray cutoff: `B = blocked_directions / directions`.
4. Retain points with `B >= min_buriedness`; group connected points using 6-neighbor adjacency.
5. By default, smooth the buriedness over grid neighbors for seed detection, select separated high-buriedness maxima, and perform a seeded, connectivity-preserving 3D watershed split of each connected region. The original buriedness values remain in the exported voxels.\n6. Rank candidates by voxel-count volume, with peak buriedness as the tie breaker.

**Caveats:** voxel volume is approximate and grid-dependent; a connected cluster may contain
multiple nearby pockets; the initial watershed segmentation is heuristic and can also over-split broad regions. Small basins below `--min-points` are currently discarded, so segmented voxel totals may be smaller than unsegmented totals. Closed internal cavities are not yet distinguished from accessible
pockets, and no biochemical druggability score is provided. These results must not be treated
as validated binding-site predictions.

## Options

| Option | Default | Meaning |
|---|---:|---|
| `--chains A,B` | all | include only selected protein chains |
| `--spacing` | 1.0 Å | grid pitch |
| `--ray-length` | 10 Å | maximum extent of rays |
| `--directions` | 42 | approximately uniform Fibonacci directions |
| `--min-buriedness` | 0.5 | minimum fraction of blocked directions |
| `--surface-distance` | 5 Å | candidate distance from nearest atomic surface |
| `--min-points` | 30 | minimum connected voxel count |
| `--top` | all | maximum number of candidate clusters exported |\n| `--no-segmentation` | off | revert to original 6-connected components |\n| `--seed-spacing` | 6 Å | minimum separation of watershed seeds |\n| `--seed-prominence` | 0.10 | minimum smoothed buriedness above the cutoff for a seed |\n| `--max-seeds` | 64 | maximum seed count in one connected region |
| `--pdb` | omitted | PDB pseudoatoms for buried voxels; B-factor = buriedness × 100 |
| `--output` | stdout | JSON summary |
| `--points` | omitted | voxel CSV for external 3D inspection |

The JSON includes the number of selected atoms, sampled voxels, clusters and, for each
exported cluster, its center, approximate volume, peak buriedness and voxel count.

## PDB visualization

Use `--pdb pockets.pdb` to export voxel centers as HETATM pseudoatoms named `C` in residue `PKT`, chain `P` (residue sequence = pocket ID). Occupancy is 1.00; the PDB temperature-factor column stores **100 × buriedness**, therefore 0–100 rather than 0–1. This is a visualization property, **not** an experimental atomic displacement factor. Color using the B-factor and display pseudoatoms as dots or spheres. The exporter rejects more than 99,999 pseudoatoms, more than 9,999 pockets, or coordinates outside standard fixed-column PDB ranges, rather than silently wrapping identifiers. JSON and CSV remain available for large exports.

## Segmentation comparison (experimental)

Compare the new default segmentation against the original clustering with:

```bash
ligint-pocket 1HSG.pdb --output pockets_split.json --pdb pockets_split.pdb
ligint-pocket 1HSG.pdb --no-segmentation --output pockets_connected.json --pdb pockets_connected.pdb
```

The separation parameters have **not** been calibrated yet against known binding sites. Do not equate the number of watershed basins with the number of experimentally established pockets.
