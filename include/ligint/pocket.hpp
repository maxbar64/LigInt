// Experimental geometric pocket search by directional buriedness. Independent of Qt.
#pragma once
#include "ligint/ligint.hpp"
#include <string>
#include <vector>

namespace ligint {
struct PocketAtom {
    V3 pos{0,0,0};
    double radius = 1.7;
    std::string chain, residue;
};
struct PocketOptions {
    double spacing = 1.0;
    double solventProbeRadius = 1.4; // radius of solvent probe, Angstrom
    double rayLength = 10.0;       // Angstrom
    double maxSurfaceDistance = 5.0;
    int morphologyRadius = 0; // 0 disabled; 1=6-neighbor cross, 2=radius-2 sphere
    std::string morphology = "none"; // none, opening, closing, open-close
    int regrowSteps = 0; // 0 disabled, 1..5 six-neighbor layers after opening/watershed
    double neckClearance = 0.0; // A; 0 disabled, otherwise cut narrow centers before clustering
    double minBuriedness = 0.5;    // fraction of blocked directions
    int directions = 42;
    std::size_t maxGridPoints = 4000000;
    std::size_t minClusterPoints = 30;
    bool segment = true;             // split large connected regions by buriedness maxima
    double seedSpacing = 6.0;        // minimum distance between distinct maxima, Angstrom
    double seedProminence = 0.10;    // minimum buriedness above the connectivity threshold
    std::size_t maxSeeds = 64;
    bool mergeBasins = true;
    double mergeSaddleDrop = 0.08;
    double mergeMaxSeedDistance = 10.0; // Angstrom; prevent long-distance chaining
    double mergeMinContactArea = 8.0; // minimum shared voxel-face area, Angstrom^2
    std::vector<std::string> chains; // empty: all; exact chain IDs
};
struct PocketPoint {
    V3 pos{0,0,0};
    float buriedness = 0;
    float solventPathDepth = -1; // shortest solvent-accessible grid path from exterior, Å; -1 if unreachable
};
struct PocketCandidate {
    V3 center{0,0,0};
    double volume = 0;             // voxel-count approximation, Angstrom^3
    double peakBuriedness = 0;
    std::size_t pointCount = 0;
    double meanBuriedness = 0;
    double highBuriednessFraction = 0; // fraction with buriedness >= 0.8
    double exposedFaceArea = 0;      // voxel boundary area, Angstrom^2
    double compactness = 0;          // 36*pi*V^2 / A^3 (voxel approximation)
    double elongation = 0;           // bounding-box major/minor extent
    double meanAtomSurfaceDistance = 0; // nearest atom center minus van der Waals radius
    double accessibilityProbe08 = 0; // reachable voxel fraction for radius 0.8 A
    double accessibilityProbe20 = 0; // reachable voxel fraction for radius 2.0 A
    double solventAccessibleFraction = 0; // fraction reachable by water-probe center
    double meanSolventPathDepth = 0;  // shortest connected path to grid exterior, Angstrom
    double maxSolventPathDepth = 0;
    double heuristicScore = 0;       // experimental ranking, not druggability
    std::vector<PocketPoint> points;
};
struct PocketSearchResult {
    std::vector<PocketCandidate> pockets;
    std::size_t sampledPoints = 0;
    std::size_t selectedAtoms = 0;
};
std::vector<PocketAtom> readPocketAtoms(const std::string& path, StructureFormat format = StructureFormat::Auto);
PocketSearchResult findPockets(const std::vector<PocketAtom>& atoms,
                               const PocketOptions& options = PocketOptions());
} // namespace ligint
