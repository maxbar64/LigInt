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

### Independent pair candidates (greedy-order diagnostic)

The grouping report also emits `pair_candidates`: every face-contacting pair passing the same contact and pair-volume thresholds, **without** excluding pairs because either member has already joined another group. This is diagnostic only, not a merge instruction. For 1HSG, the contact graph has #1–#2 with 90 shared faces and a pair volume of 1,297 Å³; it qualifies under the default 8-face, 0.02 normalized-contact and 1,800 Å³ pair-volume gates even though the greedy partition separates them. Explicitly compare pair candidates and final groups to detect order dependence. The next scientifically justified stage is a ligand-blind score for interface buriedness/clearance and independent crystallographic validation, rather than hard-coding this pair.

## Experimental C++ inter-pocket geometry and maximin throat

The independent C++17 module `include/ligint/pocket_interfaces.hpp` evaluates adjacent predicted pocket pairs without any ligand information. The CLI opt-in `--interfaces FILE` writes a separate JSON interface report with shared voxel faces, mean buriedness at the face, mean and maximum VdW clearance, and **maximin bottleneck clearance**. For each neighboring pair the widest path is traced along six-connected voxels belonging only to the two pocket candidates, between their respective maximum-clearance sites. This provides a reproducible geometric diagnostic, **not** a validated rigid-ligand passage radius, a solvent-channel width, or an automatic merging decision. It does not alter pocket generation or ranking.

```bash
ligint-pocket 1HSG.pdb --chains A,B --rank-by heuristic \\
  --morphology opening --morph-radius 1 \\
  --output opening.json --pdb opening.pdb \\
  --interfaces 1hsg_opening_interfaces.json
```

The C++ interface JSON uses pocket indices matching the same invocation's sorted CLI pocket JSON and PDB; regenerate all three files together if flags change. A second run on 1STP and 1IEP with the same settings is necessary to assess whether the interface descriptors generalize. The present calculation finds neighbors only where exported voxel centers share faces; gaps erased by opening cannot be traversed until the full pre-opening buriedness field is exposed by a later revision. The current implementation is a prototype and computationally expensive on large pair sets; no performance or validation claims are made.

## Optional overlapping pocket group beam search (prototype)

The C++17 header `include/ligint/pocket_groups.hpp` adds ligand-blind, non-destructive beam search. Each input pocket remains available; overlapping candidate sets of up to four subpockets are enumerated from measured interface adjacency, subject to a 1,800 Å³ combined-volume limit and a 0.5 Å minimum throat clearance. All original subpockets are expanded before beam pruning; subsequent layers retain up to 32 candidate groups each. The top 30 scored entries, including singleton pockets, are exported using `--groups FILE`. The ranking uses a transparent, **uncalibrated experimental** combination of normalized face contact, boundary buriedness, mean interface clearance, widest-path bottleneck, singleton heuristics and volume penalty. Groups remain candidate *unions*, not physically validated cavities. No ligand positions are used for discovery or scoring. The `--groups` option does not modify the original result.

```bash
ligint-pocket 1HSG.pdb --chains A,B --rank-by heuristic \\
  --morphology opening --morph-radius 1 \\
  --output opening.json --pdb opening.pdb \\
  --interfaces 1hsg_opening_interfaces.json \\
  --groups 1hsg_opening_groups.json
```

Pocket IDs in these JSON reports correspond to the same sorted `opening.json` invocation. Re-run all outputs together when parameters change. Validate ranking independently with 1STP and 1IEP before choosing fixed weights or enabling group ranking by default.

## Diversified Top-K site candidates and group coverage benchmark

The opt-in `--diverse-groups FILE` output evaluates up to 250 ranked candidates and selects ten via greedy maximal marginal relevance: `selection_score = raw_group_score - 0.35 * max_voxel_Jaccard(previously_selected)`. Candidate voxel sets are unions of disjoint original subpockets, so weighted Jaccard on member IDs uses their voxel volumes. This **does not modify** `--groups` output or the original pocket ranking. The overlap penalty and top-10 limit are experimental, not calibrated.

```bash
./build-pocket/src/ligint-pocket 1HSG.pdb --chains A,B \
  --rank-by heuristic --morphology opening --morph-radius 1 \
  --output opening.json --pdb opening.pdb \
  --groups 1hsg_opening_groups.json \
  --diverse-groups 1hsg_diverse_groups.json

python3 tools/pocket_group_benchmark.py \
  --crystal 1HSG.pdb --ligand MK1 --ligand-chain B --ligand-resseq 902 \
  --pockets-pdb opening.pdb --groups-json 1hsg_opening_groups.json \
  --output 1hsg_groups_benchmark.json
python3 tools/pocket_group_benchmark.py \
  --crystal 1HSG.pdb --ligand MK1 --ligand-chain B --ligand-resseq 902 \
  --pockets-pdb opening.pdb --groups-json 1hsg_diverse_groups.json \
  --output 1hsg_diverse_benchmark.json
```

The benchmark uses the ligand **only after prediction**, computes atomwise proximity to the union of group voxel sites, and records first full-coverage rank and individual group coverage. Run equivalent independent benchmarks for 1STP and 1IEP before any score calibration.

## External validation: 1STP and 1IEP

Run `python3 tools/pocket_external_benchmark.py --binary ./build-pocket/src/ligint-pocket --output-dir pocket_external_results` on Linux or Windows (provide the relevant executable). The script downloads `1STP.pdb` and `1IEP.pdb` from RCSB if absent, discovers every heavy-atom BTN or STI ligand instance directly from fixed-width PDB records, and performs *independent* protein-chain-specific runs for each instance. It retains the 1HSG baseline settings (`opening`, radius 1, heuristic rank, 3 Å post-hoc coverage) and produces per-instance raw and diversified group reports plus `pocket_external_results/external_summary.json`. Ligands are read only by the evaluation script, never passed to the predictor. Choosing a single protein chain rather than the entire biological assembly can alter accessible regions; interpret this as a chain-local benchmark and record possible assembly effects. If the PDB source format or chain mapping differs, inspect its records rather than silently selecting an arbitrary ligand. Do not tune weights on these validation structures.

## 1IEP A/B oracle coverage diagnostic

For a retrospective, ligand-informed **audit only**, run:

```bash
python3 tools/pocket_oracle_audit.py \\
 --crystal pocket_external_results/1IEP.pdb \\
 --pockets-a pocket_external_results/1IEP_A_201/opening.pdb \\
 --pockets-b pocket_external_results/1IEP_B_202/opening.pdb \\
 --output 1iep_oracle_audit.json
```

The report measures the coverage of the union of *all* original subpockets, missing ligand atom names, per-pocket contributions, and a greedy (not exact-minimum) set-cover diagnostic. This distinguishes insufficient initial voxel discovery from downstream grouping/ranking failures. Do **not** use the crystal ligand during prediction or to tune a score on the same validation structures. Note: an original subpocket larger than the 1,800 Å³ group cap is currently omitted from candidate generation entirely (observed 1IEP-A pocket #13, 2,090 Å³); treat this as a search-cap design issue, not evidence the pocket is irrelevant.

## Oversized singleton correction (1IEP-A)

Independent oracle audit showed 1IEP-A's existing singleton #13 (2,090 Å³) covers 37/37 STI heavy atoms within 3 Å, but had been removed by the 1,800 Å³ group limit. **All discovered singleton pockets are now retained**, regardless of their size. The `maxVolume` cap and volume-growth penalty still apply to *multi-pocket unions* only. The full ranking and diversified outputs continue to use the identical ligand-independent heuristic. Re-run the 1IEP external benchmark to measure whether the retained singleton reaches the final Top-K; coverage alone does not establish its ranking.

## Candidate retention audit and singleton quota (experimental)

The `--group-audit FILE` option exports all singleton pockets with their **full, untruncated geometric rank**, score, original volume, inclusion in the normal Top-30, inclusion in the 250-entry diversified pool, and final diversified position. This explicitly diagnoses search pruning without using a crystal ligand. The diversified Top-10 now reserves its first two choices for singleton pockets; the remaining eight follow the existing overlap-penalized ranking. The plain `--groups` ranking stays unchanged. This quota is a diagnostic policy, **not validated to improve Top-K coverage**; compare it against previous results before adopting as default.

```bash
ligint-pocket pocket_external_results/1IEP.pdb --chains A \\
  --rank-by heuristic --morphology opening --morph-radius 1 \\
  --output opening.json --pdb opening.pdb \\
  --groups groups.json --diverse-groups diverse.json \\
  --group-audit 1iep_a_group_audit.json
```

A pocket with `raw_rank_full > 30` was not removed at discovery; it was truncated from the report. The `in_diverse_pool` flag indicates membership of the score-ranked top 250, not actual diversified selection.

## Pocket Finder 3D workspace (experimental GUI)

The `ligint-gui` upper 3D area now has **Interactions 3D** and **Pocket Finder 3D** tabs. Open a PDB or mmCIF structure with File > Open, switch to Pocket Finder 3D, configure grid spacing, buriedness threshold, ray length, maximum surface distance, solvent probe, morphology, morphology radius and minimum cluster voxels, then select **Find pockets**. The native `findPockets` calculation runs asynchronously via the existing Qt Concurrent dependency; results are sorted by the same heuristic order as `--rank-by heuristic`. Choose **Subpockets** or **Candidate groups** and select an entry to view translucent, buriedness-tinted voxel centers surrounded by faint local receptor atoms. Drag to rotate and use the mouse wheel to zoom. No Qt 3D, OpenGL, WebEngine, RDKit or external executable is required.

This is an initial diagnostic viewer, **not yet an iso-surface / triangulated pocket boundary**. Group scores remain experimental. For reproducible auditing, the CLI still provides the richer JSON exports, interfaces, and diversified group selection; the GUI currently displays the 30 raw group candidates. Chain selection, advanced segmentation controls, saved presets and group export are planned follow-ups. If the structure is changed, stale pocket results are cleared. Use matching parameters when comparing against CLI results.

## Pocket Finder 3D controls — phase 2

The experimental GUI now exposes a **protein chain selector** populated from the opened structure; select a chain before **Find pockets** to restrict the search (or leave **All chains**). The subpocket/group list supports extended selection (Ctrl-click, Shift-click) and overlays the union of selected subpocket voxels, avoiding double drawing shared members. The lightweight software renderer also supports cool/warm palettes, opacity (10–100%), and point diameter (2–12 pixels). These are display-only controls; they do not change the geometry or ranking. The camera supports mouse drag rotation and wheel zoom. Opening a different structure clears stale results.

Limitations: voxels are rendered as depth-ordered splats rather than a true isosurface mesh; full surface triangulation, ligand overlays, advanced segmentation settings and export from the GUI are future work. The GUI is still experimental and should be validated on Linux/Windows builds.

## GUI pocket surface shell and advanced controls (phase 3)

Pocket Finder 3D offers **Voxels**, **Surface**, and **Surface + Voxels**. The surface is made from exposed faces of the existing voxel grid and is rendered by the same Qt QPainter software projection, depth sorted with alpha blending; it is deliberately **not** a smoothed marching-cubes mesh. Surface appearance depends on grid spacing; lower spacing gives a more detailed stepped envelope at additional cost. Existing palette and opacity controls apply to the shell.

The **Ligand overlay** selector reads HETATM residues from the currently opened local **PDB** file and draws their heavy-atom coordinates in yellow. It is strictly a visual diagnostic; ligands are never passed into pocket discovery, grouping or scoring. mmCIF structures continue to work for pocket discovery but do not yet have a ligand overlay in this phase.

The GUI also exposes seed separation, prominence, maximum seeds, basin merge distance and saddle drop, passed to the existing native core `PocketOptions` fields. The original CLI benchmark settings are not changed. Visual shell rendering and interactive controls are experimental; test with medium-sized structures before using large structures or small grid spacing, as projected voxel faces can be expensive.

## GUI display-only noise removal (phase 4)

The Pocket Finder 3D workspace has two independent **display-only** controls:

- `Remove isolated voxels` (default 1; 0 disables): the minimum number of occupied 6-face neighbors for a displayed voxel.
- `Min visible island` (default 12; 1 disables): hide the resulting 6-connected islands smaller than the given number of voxel centers.

Processing applies to a temporary copy of the union of currently selected pocket voxel points, before drawing either splats or exposed faces. The top overlay reports visible and hidden voxel counts. Crucially, this does not change pocket discovery, scoring, segmentation, grouping, benchmark coverage or exports. Set the controls to 0 and 1, respectively, to recover the unfiltered display. Increasing the thresholds aggressively can hide genuine narrow throats and small chemical pockets; the default values are conservative starting points rather than validated scientific thresholds.

## Pocket Finder surface relaxation (phase 5)

**Surface smoothing** (0–4, default 2) applies a display-only, shared-corner occupancy-centroid relaxation to the exposed voxel faces. Zero reproduces the previous sharp voxel shell. Increasing the value reduces blocky local corners; the relaxation uses the six-connected island and neighbor gates configured above and does not mutate `PocketCandidate::points`. It does not implement true Marching Cubes, curvature flow, or a density isosurface; some faceting may remain, especially on sparse or narrow pocket regions. Benchmark results and crystal ligand coverage remain unchanged. Compare levels 0/1/2/3/4 on 1IEP-A pocket 13 before selecting a GUI default for routine use. Rendering geometry is generated inside paintEvent, so large cavities may incur increased redraw time while rotating.

## Smooth triangulated isosurface (experimental)

The 3D rendering menu includes **Smooth isosurface** alongside **Voxels**, **Voxel shell**, and **Shell + Voxels**. This new option rasterizes an actual triangular isosurface generated from a local voxel occupancy scalar field. A separable-binomial-style 3×3×3 weighted Gaussian blur is applied from 0 to 4 passes (the existing **Surface smoothing** control); the surface is extracted using marching tetrahedra. The pre-existing display-only neighbor/island filters apply before scalar-field construction. The result is triangulated, distinct from deforming the old blocky shell's corners. The mesh is cached across camera rotation and zoom, and invalidated on parameter/selection/new-search changes. There is a 3.5-million-cell safety cap for large selections; beyond it, the smooth mesh is skipped. All transformations operate on display-only copies: pocket scoring and docking inputs are not modified. Because mesh colors are currently unified per selection and rasterized with painter sorting, this is a first experimental geometry path rather than final high-end rendering. Compare **Voxel shell** and **Smooth isosurface**, with smoothing 1–3, on 1IEP-A pocket #13. Confirm Qt GUI compilation and test on Linux and Windows.

## Pocket 3D lighting (experimental)

The **Smooth isosurface** render mode now shades triangles with camera-space diffuse lighting and an adjustable specular highlight. The **Diffuse light** control ranges from 0–100% (default 75%); **Specular highlight** ranges from 0–100% (default 30%). Triangle normals are computed as normalized cross products of the cached mesh vertices, with two-sided diffuse response to reduce black patches from locally inconsistent marching-tetrahedra winding. A fixed camera-space light direction updates each frame on rotation without recomputing the mesh. This is a preliminary **flat-shaded** treatment: it is not per-vertex smooth normals, physically based shading, true per-pixel normal interpolation, or order-independent transparency. Colors, geometry and light controls are independent of scientific pocket calculations. Test on the 1IEP-A pocket #13 and compare diffuse 0/75/100 and highlight 0/30/80; if triangle facets remain conspicuous, the subsequent rendering task should be mesh vertex-normal accumulation and improved depth-aware transparency.
