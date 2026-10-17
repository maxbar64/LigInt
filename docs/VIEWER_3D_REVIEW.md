# 3D Viewer Review

## Scope

Reviewed the current `main` snapshot at commit `69f626340f0dedc7ecd957580bef989d1f984ca1`, focusing on `gui/ViewerWidget.{h,cpp}`, its integration in `gui/MainWindow.cpp`, the GUI smoke test, and `docs/ARCHITECTURE.md`.

## Current design

- `ViewerWidget` renders ligand atoms/bonds, nearby contact atoms/bonds, interaction vectors, labels, distances and a legend with Qt `QPainter`.
- The camera is an orthographic projection. The initial orientation uses the principal axes of ligand and interaction coordinates; mouse drag rotates, right/middle drag pans, the wheel zooms, and double click resets.
- Atom and bond primitives are depth-sorted and use depth-dependent opacity. The drawing code is self-contained in the GUI and does not require OpenGL.
- `MainWindow` exposes labels, residues, distances and legend toggles. Selecting a row in the interaction table highlights that interaction in the viewer.
- `renderImage()` shares the screen renderer and supports image export and offscreen testing.

This is a compact, portable baseline for an interaction-focused view. The current renderer is a 2D projection of a 3D scene; it is not a GPU-backed molecular graphics engine.

## Review findings

1. **Interactions are painted after the depth-sorted molecular scene.** Interaction lines, water markers and ring outlines are drawn after all atoms and bonds, without depth testing. A line behind a foreground atom can therefore remain visible through it. This is a visual occlusion limitation; confirm whether that is intentional before changing the visual convention.
2. **The initial framing is interaction-centred and deliberately bounded.** The camera bounds use ligand atoms and interaction endpoints, then cap the contribution of contact atoms. This keeps the ligand prominent, but it will not support a complete pocket or whole-protein context view without a separate framing mode.
3. **The viewport does not pick objects.** Highlighting is driven by the table selection; atoms and interaction marks cannot currently be selected directly in the view.
4. **The view controls are global toggles only.** There is no projection selector, camera preset, atom-name display, residue-selection control, or independent visibility for interaction lines and molecular components.
5. **Rendering work is repeated on every paint.** Coordinates are projected multiple times and the depth-sorted item list is rebuilt every frame. This is reasonable for the current local scene, but should be measured before expanding to larger residue or protein selections.
6. **GUI rendering coverage is smoke-level.** The test checks successful rendering and image diversity, but does not assert projection geometry, occlusion, picking, camera state, or visibility combinations.

No correctness regression was demonstrated by this review. Findings 1 and 2 are current design constraints that should guide the feature design and tests.

## Branch implementation

The first branch increment adds a dedicated **Graphics options** panel next to the viewer. It groups the existing residue, label, distance and legend toggles under **Structure** and **Interactions**, and puts **Reset view** under **Camera**. The previous toolbar controls have been moved into this panel. The GUI smoke test checks that the panel is visible, the controls start in the expected state, and a toggle can be changed. The manual source and generated HTML describe the new location.

## Suggested implementation sequence

1. Agree on the first user-visible capability and its default behavior.
2. Add a small camera/scene state API to `ViewerWidget`, preserving the present orthographic view and image export.
3. Implement the selected capability with focused GUI tests for geometry, visibility and reset behavior.
4. Extend the smoke test with image-level assertions only where stable; prefer state and projection assertions over fragile pixel snapshots.
5. Revisit caching or a GPU renderer only if profiling shows the expanded scene needs it.

Potential first increments are perspective projection, viewport selection of interactions/atoms, or a selectable pocket/residue context. These have different interaction and rendering consequences, so this branch leaves that choice open.
