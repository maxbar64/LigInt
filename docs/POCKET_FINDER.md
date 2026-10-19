# Experimental buriedness pocket finder

`ligint-pocket` is a standalone CLI for ligand-independent geometric pocket discovery.
It is built with `LIGINT_BUILD_CLI=ON`. It does **not** modify the regular `ligint` interaction analysis.

```bash
ligint-pocket protein.pdb --chains A,B --output pockets.json --points voxels.csv
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
5. Rank candidates by voxel-count volume, with peak buriedness as the tie breaker.

**Caveats:** voxel volume is approximate and grid-dependent; a connected cluster may contain
multiple nearby pockets. Closed internal cavities are not yet distinguished from accessible
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
| `--top` | all | maximum number of candidate clusters exported |
| `--output` | stdout | JSON summary |
| `--points` | omitted | voxel CSV for external 3D inspection |

The JSON includes the number of selected atoms, sampled voxels, clusters and, for each
exported cluster, its center, approximate volume, peak buriedness and voxel count.
