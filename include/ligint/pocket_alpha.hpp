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
    std::size_t minSpheres=3, maxCandidates=200000;
};
struct AlphaSphere { V3 center{0,0,0}; double radius=0; };
struct AlphaPocket { V3 center{0,0,0}; double volume=0; std::vector<AlphaSphere> spheres; };
std::vector<AlphaPocket> findAlphaPockets(const std::vector<PocketAtom>& atoms,
                                          const AlphaPocketOptions& options={});
}
