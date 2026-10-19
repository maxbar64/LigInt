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

## Optional 3D morphology and narrow-connection pruning

Experimental post-buriedness processing operates **before** connected components
and watershed. Both are disabled by default to preserve the reference results.

```bash
ligint-pocket 1HSG.pdb --chains A,B --morphology opening --morph-radius 1 \
  --output 1hsg_opening.json --pdb 1hsg_opening.pdb
ligint-pocket 1HSG.pdb --chains A,B --neck-clearance 0.8 \
  --output 1hsg_neck.json --pdb 1hsg_neck.pdb
```

Options:

- `--morphology none|opening|closing|open-close`: opening erodes then dilates the selected voxel mask; closing dilates then erodes; open-close applies both, using a spherical voxel neighborhood.
- `--morph-radius 0|1|2`: discrete structuring-element radius in **grid voxels**. Default 0 disables morphology, including when a non-none mode is selected.
- `--neck-clearance A`: remove voxel centers with nearest atom-surface clearance less than A Å (default 0, disabled), without editing the original buriedness measurements.

The original scalar buriedness is not smoothed or overwritten. A closing operation may recover sampled, previously below-threshold voxel centers, using their **measured** buriedness; never assign new values to unsampled positions or to points inside atomic volume. Consequently, a recovered voxel can have a buriedness below `--min-buriedness`. Morphological opening and neck removal can reduce volume and ligand coverage; compare them against the baseline rather than assuming an improvement.

The neck-clearance filter is a **hard clearance mask**, not yet a medial-axis or distance-transform bottleneck segmentation; its effect is to remove narrow space globally. Spherical elements at 1 Å grid resolution are strongly discretized. Grid-pitch sensitivity and preservation of crystal-ligand coverage are required before enabling either filter by default.

Suggested controlled comparison: baseline, opening radius 1, and clearance thresholds 0.5/0.8 Å on **both chains** of 1HSG, plus 1STP and 1IEP. Record coverage within 2 and 3 Å, retained voxels, the number of pockets, and the largest region's fraction of total volume.

## Bounded regrowth from opening cores (experimental)

Use `--morphology opening --morph-radius 1 --regrow-steps 2` to restore at most two 6-connected grid layers of measured, above-threshold buriedness from the **original unfiltered score** after segmentation and basin merging. All surviving final basins expand together. Equidistant voxels contested by multiple basins are deliberately left unassigned; no newly measured or synthetic buriedness is generated. Pocket centers, volumes and descriptors are recomputed after regrowth. The baseline remains `--regrow-steps 0` (disabled). This first version does not recover below-threshold or unsampled positions, and the 6-neighbor growth layers are grid-dependent rather than true Euclidean distances. `--regrow-steps` requires opening or open-close with nonzero morph radius, and is capped at five layers.

Suggested 1HSG A,B comparison (repeat with 1STP and 1IEP):

```bash
ligint-pocket 1HSG.pdb --chains A,B --rank-by heuristic --morphology opening --morph-radius 1 --regrow-steps 2 --output 1hsg_regrow2.json --pdb 1hsg_regrow2.pdb
python3 tools/pocket_benchmark.py --crystal 1HSG.pdb --ligand MK1 --ligand-chain B --ligand-resseq 902 --pockets-json 1hsg_regrow2.json --pockets-pdb 1hsg_regrow2.pdb --output 1hsg_regrow2_benchmark.json
```

Check ligand coverage, largest-pocket volume, Top-1/3/5 ranking and whether regrowth makes previously small, isolated pockets too broad. Do not make regrowth the default from the 1HSG result alone.

## Crystal-ligand atom and pocket-contact audit

The independent evaluator `tools/pocket_contact_audit.py` reports the closest pocket to each heavy ligand atom, all pockets within 3 Å and a **ligand-blind** six-connected contact graph between distinct exported pockets. The graph uses a one-voxel face adjacency on the PDB pseudoatom coordinates. The ligand is used only for the audit after pocket discovery, never to produce or merge pockets.

```bash
python3 tools/pocket_contact_audit.py --crystal 1HSG.pdb --ligand MK1 \\
  --ligand-chain B --ligand-resseq 902 \\
  --pockets-pdb opening.pdb --spacing 1.0 \\
  --output 1hsg_opening_contacts.json
```

On the 1HSG opening reference, four of 45 MK1 heavy atoms are outside the best single pocket's 3 Å coverage: C4, C5, C6, O2. All four lie close to *other* pockets. Contact alone is NOT sufficient grounds for merging two basins: a future ligand-independent merge criterion must also check saddle/buriedness continuity, neck width, and excessive global component sizes, and be validated on 1STP and 1IEP.

## Ligand-blind hierarchical candidate groups (experimental)

The stand-alone program `tools/pocket_group_candidates.py` reads an existing pseudoatom PDB and proposes **non-destructive** groups of face-adjacent pocket subregions. The original predicted regions, geometry and ranking remain unchanged. Edges are ordered by shared-face fraction `faces/(6*min(voxels_A, voxels_B))`; greedy grouping is constrained by a maximum group voxel-volume and a maximum number of subregions. All grouping parameters are experimental and need multi-crystal validation. This prototype does **not yet** implement saddle/buriedness continuity or a throat-radius check, and therefore is not a final merge classifier.

```bash
python3 tools/pocket_group_candidates.py --pockets-pdb opening.pdb \\
  --spacing 1 --min-contact-faces 8 --min-contact-fraction 0.02 \\
  --max-group-volume 1800 --max-group-members 4 \\
  --output 1hsg_group_candidates.json
python3 -m unittest discover -s tools -p 'test_pocket_group_candidates.py'
```

In the 1HSG opening reference (55 subpockets), initial contact+volume grouping produces ten multiregion candidate groups. One such group contains pockets 2, 3, 14 and 43 (about 1,374 Å³); it **does not** automatically select the crystal-evaluation pair 1+2. This negative result is valuable: do not tune the grouping to MK1 or interpret a contact heuristic as proof of a binding site. Validate ligand coverage and pocket precision post-hoc on independent 1STP and 1IEP complexes before modifying the core predictor.
