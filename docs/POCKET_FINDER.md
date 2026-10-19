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
5. By default, smooth the buriedness over grid neighbors for seed detection, select separated high-buriedness maxima, and perform a seeded, connectivity-preserving 3D watershed split of each connected region. The original buriedness values remain in the exported voxels.\n6. Before exporting, merge adjacent watershed basins whose highest shared buriedness saddle is within `--merge-saddle-drop` of the weaker peak.\n7. Rank candidates by voxel-count volume, with peak buriedness as the tie breaker.

**Caveats:** voxel volume is approximate and grid-dependent; a connected cluster may contain
multiple nearby pockets; the initial watershed segmentation is heuristic and can also over-split broad regions. Small basins below `--min-points` are currently discarded, so segmented voxel totals may be smaller than unsegmented totals. Closed internal cavities are not yet distinguished from accessible
pockets, and no biochemical druggability score is provided. These results must not be treated
as validated binding-site predictions.

## Options

| Option | Default | Meaning |
|---|---:|---|
| `--chains A,B` | all | include only selected protein chains |
| `--spacing` | 1.0 Å | grid pitch |\n| `--solvent-probe` | 1.4 Å | probe-center excluded volume radius |
| `--ray-length` | 10 Å | maximum extent of rays |
| `--directions` | 42 | approximately uniform Fibonacci directions |
| `--min-buriedness` | 0.5 | minimum fraction of blocked directions |
| `--surface-distance` | 5 Å | candidate distance from nearest atomic surface |
| `--min-points` | 30 | minimum connected voxel count |
| `--top` | all | maximum number of candidate clusters exported |
| `--rank-by` | volume | `volume` (legacy) or `heuristic` (experimental score) |\n| `--no-segmentation` | off | revert to original 6-connected components |\n| `--seed-spacing` | 6 Å | minimum separation of watershed seeds |\n| `--seed-prominence` | 0.10 | minimum smoothed buriedness above the cutoff for a seed |\n| `--max-seeds` | 64 | maximum seed count in one connected region |\n| `--merge-saddle-drop` | 0.08 | maximum buriedness drop from weaker peak to inter-basin saddle for merging |\n| `--merge-max-distance` | 10 Å | maximum span between any two member seeds in a merged basin |\n| `--merge-min-contact` | 8 Å² | minimum shared voxel-face area for an eligible merge |\n| `--no-basin-merge` | off | preserve all watershed basins |
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

The merging stage uses adjacency of 6-connected voxels and the strongest shared buriedness saddle. Results are heuristic; vary `--merge-saddle-drop` (e.g. 0.04, 0.08, 0.12) and compare pocket sizes and experimental ligand coverage. The merged candidate's reported buriedness and original voxel values are preserved.

Geometric merge gates: each candidate pair must share at least `--merge-min-contact` square angstroms of voxel faces, and all maxima in the proposed merged group must lie within `--merge-max-distance` of one another. These restrictions are designed to prevent chaining across the entire protein; the contact-area estimate is grid dependent.

## Geometric characterization (experimental)

The JSON summary now includes for each retained region:

- `mean_buriedness`: mean of original directional buriedness samples (0 to 1).
- `high_buriedness_fraction`: fraction of points with buriedness >= 0.8.
- `boundary_area_a2`: exposed 6-neighbor voxel-face area; grid dependent.
- `compactness`: dimensionless isoperimetric quotient `36*pi*volume^2/area^3`; computed from the voxel area, grid dependent.
- `elongation`: longest/shortest axis-aligned bounding-box extent. Depends on structure orientation.
- `mean_surface_distance_a`: mean nearest atomic van der Waals surface distance, clipped at zero. This is **not** solvent-accessible depth and should not be presented as such.
- `heuristic_score`: **unvalidated** composite using volume saturation, mean buriedness, highly buried fraction and compactness. Not an affinity, druggability or probability estimate.

The default pocket IDs remain sorted by **volume**, as in earlier revisions. Use `--rank-by heuristic` to sort JSON, CSV and PDB pocket IDs by decreasing `heuristic_score`. No ligand coordinates or experimentally annotated binding sites enter the score.

To compare 1HSG with the earlier segmentation, run the same parameters on chain A, export JSON and compare `pockets[].heuristic_score` and geometric descriptors while keeping the original PDB/CSV voxel geometry unchanged.

## Solvent access and path-length proxies

After pocket segmentation, a probe-center exclusion grid is computed by inflating each atom's van der Waals radius by `--solvent-probe` (default 1.4 Å). A 6-neighbor flood fill from the grid boundary identifies probe-center-accessible grid sites. For each pocket, JSON includes:

- `solvent_accessible_fraction`: fraction of pocket voxel centers reachable by the probe center; voxels within inflated atoms count as inaccessible.
- `mean_solvent_path_depth_a`: average shortest **6-neighbor grid path to the exterior boundary** for accessible pocket voxel centers, in Å; zero when none are reachable.
- `max_solvent_path_depth_a`: maximum such distance, in Å.

These are grid-based accessibility/path proxies, **not** accurate pocket-mouth depth measures. They depend on grid size/pitch, probe radius, and whether a probe fits through narrow corridors. Pockets with no reachable center have accessible fraction zero and zero reported depth; interpret the pair jointly. The pre-existing heuristic score and volume-based default order are unchanged until benchmark validation.

## Multi-probe accessibility (experimental)

Three independent exterior-connected probe-center maps are now evaluated: 0.8 Å, the configurable `--solvent-probe` (default 1.4 Å), and 2.0 Å. The JSON includes `accessibility_probe_0_8`, `solvent_accessible_fraction`, and `accessibility_probe_2_0`. Fractions for larger probes cannot exceed those of smaller probes on the same grid (assuming the middle radius lies between the other two). A large difference is evidence of narrow accessible corridors, not by itself proof of a ligand-passable entrance. The pre-existing ranking is unchanged.

## Multi-crystal ligand-independent benchmark

Start with these experimentally determined RCSB entries, covering different binding environments:

| PDB | Protein and ligand | Ligand residue name | Why useful |
|---|---|---|---|
| 1HSG | HIV-1 protease / L-735,524 | MK1 | Dimer-interface reference; include **both A and B protein chains** |
| 1STP | Streptavidin / biotin | BTN | Deep high-affinity binding pocket |
| 1IEP | ABL kinase / imatinib | STI | Large kinase inhibitor site |

Retrieve structures from `https://files.rcsb.org/download/1HSG.pdb`, `1STP.pdb`, and `1IEP.pdb`. Confirm the desired ligand instance (auth chain and residue number) in each local PDB; the benchmark script deliberately refuses ambiguous ligand selections. Downloaded data are not stored in this repository.

Run the pocket finder **without any ligand coordinates** (it loads protein ATOM records only), then score the resulting pockets against the crystal ligand **after prediction**. For 1HSG:

```bash
ligint-pocket 1HSG.pdb --chains A,B --rank-by heuristic \\
  --output 1hsg_ab.json --pdb 1hsg_ab_pockets.pdb
python3 tools/pocket_benchmark.py --crystal 1HSG.pdb --ligand MK1 \\
  --ligand-chain B --ligand-resseq 902 \\
  --pockets-pdb 1hsg_ab_pockets.pdb --pockets-json 1hsg_ab.json \\
  --output 1hsg_ab_benchmark.json
```

For another complex, replace the crystal file and ligand code; specify its actual crystallographic chain and residue number. Use Top-1, Top-3, Top-5 and the best-pocket per-ligand heavy-atom coverage within 3 Å. **These metrics are voxel-proximity metrics, not a validated pocket-overlap or docking-success measure.** Record chain selection, crystallographic assembly, probe radii, grid parameters and ligand-instance identity for reproducibility. Do not tune ranking coefficients to the single 1HSG example.
