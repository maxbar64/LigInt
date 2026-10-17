// ligint - geometry for the 3D views: protein backbone, secondary structure, cartoon and molecular surfaces.
//
// Everything here is plain C++ (no Qt, no OpenGL): the library produces triangle meshes and the front-ends draw them.
// The backbone and the atoms are filled in Result when Options::protein3d is set (they are not needed by the command line).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "ligint/ligint.hpp"

namespace ligint {

// BackboneResidue and ProteinAtom (the protein description for the graphics) are in ligint.hpp: Result carries them.

// Secondary structure by the method of Kabsch and Sander (DSSP): hydrogen-bond energies between backbone atoms,
// n-turns for helices (3-10, alpha and pi are all 'H'), bridges and ladders for strands. Fills BackboneResidue::ss.
void assignSecondaryStructure(std::vector<BackboneResidue>& residues);

// ---------------------------------------------------------------------------------------------
// Triangle meshes
// ---------------------------------------------------------------------------------------------
using V3f = std::array<float, 3>;
using Tri = std::array<int, 3>;

struct Mesh {
    std::vector<V3f> pos;       // vertices
    std::vector<V3f> nrm;       // unit normals, pointing outwards
    std::vector<float> tag;     // cartoon: fractional index in the backbone (residue + position along the segment);
                                // surface: index of the nearest atom
    std::vector<Tri> tri;       // counter-clockwise when seen from outside
    bool empty() const { return tri.empty(); }
    void append(const Mesh& other);
};

struct CartoonOptions {
    double coilRadius = 0.28;           // loops: round tube
    double helixHalfWidth = 1.25;       // helices: flat ribbon
    double helixHalfThickness = 0.22;
    double strandHalfWidth = 0.95;      // strands: flat ribbon with an arrow head
    double strandHalfThickness = 0.20;
    double arrowHalfWidth = 1.75;
    bool tubeOnly = false;              // every residue as a tube of tubeRadius (no secondary-structure shapes)
    double tubeRadius = 0.42;
    int quality = 2;                    // 0 fast (large proteins), 1 medium, 2 high; reduced automatically for huge inputs
};

// Cartoon of the whole protein. Mesh::tag is the fractional residue index, so a front-end can color by residue,
// by secondary structure, by chain or with a smooth gradient along the chain.
Mesh buildCartoon(const std::vector<BackboneResidue>& residues, const CartoonOptions& options = CartoonOptions());

enum class SurfaceKind { SolventExcluded, SolventAccessible };

struct SurfaceOptions {
    SurfaceKind kind = SurfaceKind::SolventExcluded;
    double probe = 1.4;                 // solvent radius (Angstrom)
    double spacing = 0.5;               // grid step; raised automatically for very large systems
    int smoothing = 2;                  // iterations of Taubin smoothing of the vertices (0 = none)
};

// Molecular surface of a set of atoms; the cavities that the solvent cannot reach are removed.
// Mesh::tag is the index of the nearest atom of `atoms`.
Mesh buildSurface(const std::vector<ProteinAtom>& atoms, const SurfaceOptions& options = SurfaceOptions());

// Keeps the triangles whose center is within `radius` of one of the points (e.g. the pocket around a ligand).
Mesh clipMeshNear(const Mesh& mesh, const std::vector<V3>& points, double radius);

double meshArea(const Mesh& mesh);

}  // namespace ligint
