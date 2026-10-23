// Experimental native empty-sphere pocket candidates, independent of directional buriedness.
// Approximate Voronoi vertices generated from local atom quadruplets; not a reimplementation of fpocket.
#pragma once
#include "ligint/pocket.hpp"
#include <vector>
namespace ligint {
struct AlphaPocketOptions {
    double minRadius=2.0, maxRadius=5.5;
    double neighborCutoff=10.0;
    double clusterDistance=3.0;
    double maxClusterRadius=9.0;
    std::size_t minContactAtoms=4;
    double macroLinkDistance=5.0; // nearest sphere-center link, Å
    double macroMaxDiameter=18.0; // cap on merged microcluster center diameter, Å
    double macroMinClearance=1.0; // minimum VdW clearance on sphere-center connection, Å
    bool mergeMicroclusters=true;
    std::size_t minSpheres=3, maxCandidates=200000;
};
struct AlphaSphere { V3 center{0,0,0}; double radius=0; };
struct AlphaPocket { V3 center{0,0,0}; double volume=0; std::vector<AlphaSphere> spheres; };
struct AlphaMergeDecision {
    std::size_t first=0, second=0; // 1-based microcluster indices
    double linkDistance=0, minClearance=0, endpointClearance=0, mergedDiameter=0;
    bool accepted=false;
    const char* reason=""; // fixed static diagnostic string
};
struct AlphaMergeAudit {
    std::size_t microclusters=0;
    std::vector<AlphaMergeDecision> decisions;
};
std::vector<AlphaPocket> findAlphaPockets(const std::vector<PocketAtom>& atoms,
                                          const AlphaPocketOptions& options={},
                                          AlphaMergeAudit* audit=nullptr);
}
