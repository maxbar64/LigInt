// Experimental, ligand-independent pocket interface characterization.
#pragma once
#include "ligint/pocket.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace ligint {
struct PocketInterface {
    std::size_t first=0, second=0; // zero-based indices in supplied pockets
    std::size_t contactFaces=0;
    double contactArea=0;
    double meanBoundaryBuriedness=0;
    double meanBoundaryClearance=0;
    double maxBoundaryClearance=0;
    double bottleneckClearance=0; // widest path between the maximum-clearance grid sites
};
inline std::vector<PocketInterface> analyzePocketInterfaces(
    const std::vector<PocketCandidate>& pockets,
    const std::vector<PocketAtom>& atoms,double spacing=1.0) {
    if(!(spacing>0) || !std::isfinite(spacing))
        throw std::invalid_argument("Pocket interface spacing must be positive");
    using Key=std::array<long long,3>;
    struct Site {std::size_t pocket; float buriedness; double clearance;};
    std::map<Key,Site> sites;
    if(pockets.empty())return {};
    std::array<double,3> origin=pockets.front().points.empty()
        ? std::array<double,3>{0,0,0} : pockets.front().points.front().pos;
    for(const auto& p:pockets)for(const auto& pt:p.points)
        for(int d=0;d<3;d++)origin[d]=std::min(origin[d],pt.pos[d]);
    auto keyOf=[&](const V3& p){
        Key key{};
        for(int d=0;d<3;d++)key[d]=static_cast<long long>(std::llround((p[d]-origin[d])/spacing));
        return key;
    };
    for(std::size_t i=0;i<pockets.size();i++)for(const auto& pt:pockets[i].points){
        double clearance=std::numeric_limits<double>::infinity();
        for(const auto& atom:atoms){
            double d2=0;
            for(int k=0;k<3;k++)d2+=(pt.pos[k]-atom.pos[k])*(pt.pos[k]-atom.pos[k]);
            clearance=std::min(clearance,std::sqrt(d2)-atom.radius);
        }
        if(atoms.empty())clearance=0;
        auto inserted=sites.emplace(keyOf(pt.pos),Site{i,pt.buriedness,clearance});
        if(!inserted.second)
            throw std::invalid_argument("Duplicate voxel in pocket input");
    }
    const Key offsets[3]={Key{1,0,0},Key{0,1,0},Key{0,0,1}};
    struct Accum {std::size_t faces=0;double sumB=0,sumC=0,maxC=0;};
    std::map<std::pair<std::size_t,std::size_t>,Accum> accum;
    for(const auto& site:sites)for(const auto& direction:offsets){
        Key adjacent=site.first;
        for(int k=0;k<3;k++)adjacent[k]+=direction[k];
        auto it=sites.find(adjacent);
        if(it==sites.end() || site.second.pocket==it->second.pocket)continue;
        auto pair=std::minmax(site.second.pocket,it->second.pocket);
        auto& a=accum[{pair.first,pair.second}];
        a.faces++;
        a.sumB+=0.5*(site.second.buriedness+it->second.buriedness);
        const double c=std::min(site.second.clearance,it->second.clearance);
        a.sumC+=c;
        a.maxC=std::max(a.maxC,c);
    }
    std::vector<PocketInterface> results;
    for(const auto& edge:accum){
        const auto a=edge.first.first,b=edge.first.second;
        const auto& ac=edge.second;
        PocketInterface result;
        result.first=a;result.second=b;
        result.contactFaces=ac.faces;
        result.contactArea=ac.faces*spacing*spacing;
        result.meanBoundaryBuriedness=ac.sumB/ac.faces;
        result.meanBoundaryClearance=ac.sumC/ac.faces;
        result.maxBoundaryClearance=ac.maxC;
        // Maximin (widest) path restricted to the union of the two subpockets.
        // Start at the maximum-clearance site in A, end at the same in B.
        using Node=std::pair<double,Key>;
        std::priority_queue<Node> todo;
        std::map<Key,double> best;
        Key start{},target{};
        double startC=-std::numeric_limits<double>::infinity();
        double targetC=-std::numeric_limits<double>::infinity();
        for(const auto& cell:sites) {
            if(cell.second.pocket==a && cell.second.clearance>startC) {
                start=cell.first;startC=cell.second.clearance;
            }
            if(cell.second.pocket==b && cell.second.clearance>targetC) {
                target=cell.first;targetC=cell.second.clearance;
            }
        }
        best[start]=startC;
        todo.push({startC,start});
        const Key directions[6]={Key{1,0,0},Key{-1,0,0},Key{0,1,0},
                                 Key{0,-1,0},Key{0,0,1},Key{0,0,-1}};
        while(!todo.empty()){
            auto current=todo.top();todo.pop();
            const auto found=best.find(current.second);
            if(found==best.end() || current.first<found->second)continue;
            if(current.second==target){result.bottleneckClearance=current.first;break;}
            for(const auto& delta:directions){
                Key k=current.second;
                for(int j=0;j<3;j++)k[j]+=delta[j];
                const auto it=sites.find(k);
                if(it==sites.end() || (it->second.pocket!=a && it->second.pocket!=b))continue;
                const double quality=std::min(current.first,it->second.clearance);
                auto old=best.find(k);
                if(old==best.end() || quality>old->second){
                    best[k]=quality;todo.push({quality,k});
                }
            }
        }
        results.push_back(result);
    }
    return results;
}
} // namespace ligint
