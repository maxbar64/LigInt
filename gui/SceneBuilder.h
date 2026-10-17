// Builds the triangles and spheres of the 3D view from the analysis result and the chosen style.
// The functions only read their arguments and return new objects: the heavy ones can run on a worker thread.
#pragma once

#include <QColor>
#include <cstdint>
#include <memory>
#include <vector>

#include "Raster.h"
#include "ViewStyle.h"
#include "ligint/ligint.hpp"
#include "ligint/scene3d.hpp"

namespace gui {

struct MeshPiece {
    ligint::Mesh mesh;
    std::vector<std::uint32_t> colors;     // 0x00RRGGBB per vertex
};

struct AtomsPiece {
    ligint::Mesh sticks;
    std::vector<std::uint32_t> stickColors;
    std::vector<raster::Sphere> spheres;
};

std::uint32_t packColor(const QColor& c);
QColor backgroundColor(Background b);

// Cartoon / tube / trace of the whole protein, colored by the chosen scheme.
std::shared_ptr<const MeshPiece> makeCartoonPiece(const ligint::Result& result, const ViewStyle& style);

// Molecular surface of the protein (the whole protein or the pocket around the ligand). Slow for large proteins.
std::shared_ptr<const MeshPiece> makeSurfacePiece(const ligint::Result& result, const ligint::LigandResult& ligand, const ViewStyle& style);

// The ligand, the residues in contact and the waters of the water bridges.
std::shared_ptr<const AtomsPiece> makeLigandPiece(const ligint::LigandResult& ligand, const ViewStyle& style);
std::shared_ptr<const AtomsPiece> makeResiduePiece(const ligint::LigandResult& ligand, const ViewStyle& style);
std::shared_ptr<const AtomsPiece> makeWaterPiece(const ligint::LigandResult& ligand, const ViewStyle& style);

}  // namespace gui
