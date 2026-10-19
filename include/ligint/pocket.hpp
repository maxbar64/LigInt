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
    double spacing = 1.0;          // Angstrom
    double rayLength = 10.0;       // Angstrom
    double maxSurfaceDistance = 5.0;
    double minBuriedness = 0.5;    // fraction of blocked directions
    int directions = 42;
    std::size_t maxGridPoints = 4000000;
    std::size_t minClusterPoints = 30;
    bool segment = true;             // split large connected regions by buriedness maxima
    double seedSpacing = 6.0;        // minimum distance between distinct maxima, Angstrom
    double seedProminence = 0.10;    // minimum buriedness above the connectivity threshold
    std::size_t maxSeeds = 64;
    bool mergeBasins = true;
    double mergeSaddleDrop = 0.08; // merge when saddle is within this of lower peak       // maximum basin seeds per connected region
    std::vector<std::string> chains; // empty: all; exact chain IDs
};
struct PocketPoint {
    V3 pos{0,0,0};
    float buriedness = 0;
};
struct PocketCandidate {
    V3 center{0,0,0};
    double volume = 0;             // voxel-count approximation, Angstrom^3
    double peakBuriedness = 0;
    std::size_t pointCount = 0;
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
