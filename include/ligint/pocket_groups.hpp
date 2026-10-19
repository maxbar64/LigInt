// Experimental overlap-permitting beam search over ligand-blind pocket interfaces.
#pragma once
#include "ligint/pocket_interfaces.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <numeric>
#include <stdexcept>
#include <vector>
namespace ligint {
struct PocketGroup {
    std::vector<std::size_t> members; // zero-based indices
    double volume=0, score=0, weakestThroat=0;
};
inline std::vector<PocketGroup> rankPocketGroups(
 const std::vector<PocketCandidate>& pockets,
 const std::vector<PocketInterface>& edges,
 std::size_t beamWidth=32,std::size_t topK=30,std::size_t maxMembers=4,
 double maxVolume=1800,double minThroat=0.5){
    if(beamWidth==0||maxMembers==0||!(maxVolume>0)||!std::isfinite(maxVolume)
       ||!std::isfinite(minThroat))
        throw std::invalid_argument("Invalid group search options");
    const auto n=pockets.size();
    struct E {std::size_t other;double quality,throat;};
    std::vector<std::vector<E>> adjacency(n);
    for(const auto& e:edges){
        if(e.first>=n||e.second>=n||e.first==e.second)
            throw std::invalid_argument("Invalid pocket interface indices");
        if(e.bottleneckClearance<minThroat || e.contactFaces==0)continue;
        // Normalization against the smaller pocket makes this scale independent.
        const double fraction=static_cast<double>(e.contactFaces)/
            (6.0*std::max<std::size_t>(1,std::min(pockets[e.first].pointCount,pockets[e.second].pointCount)));
        const double quality=0.30*std::min(1.0,fraction/0.08)
           +0.25*std::clamp(e.meanBoundaryBuriedness,0.0,1.0)
           +0.30*std::clamp(e.bottleneckClearance/2.5,0.0,1.0)
           +0.15*std::clamp(e.meanBoundaryClearance/2.0,0.0,1.0);
        adjacency[e.first].push_back({e.second,quality,e.bottleneckClearance});
        adjacency[e.second].push_back({e.first,quality,e.bottleneckClearance});
    }
    auto score=[&](PocketGroup& g){
        double total=0,edgeSum=0;std::size_t contacts=0;
        for(auto id:g.members)total+=pockets[id].volume;
        g.volume=total;
        g.weakestThroat=0;
        // Best possible spanning tree edges in the induced candidate graph.
        std::set<std::size_t> reached{g.members.front()};
        double minimum=1e9;
        while(reached.size()<g.members.size()){
            double best=-1;std::size_t next=n;double throat=0;
            for(auto a:reached)for(const auto& e:adjacency[a])
                if(std::binary_search(g.members.begin(),g.members.end(),e.other)
                   && !reached.count(e.other) && e.quality>best){
                    best=e.quality;next=e.other;throat=e.throat;
                }
            if(next==n){g.score=-1e9;return;}
            reached.insert(next);edgeSum+=best;++contacts;
            minimum=std::min(minimum,throat);
        }
        g.weakestThroat=contacts?minimum:0;
        const double pocketPrior=std::accumulate(g.members.begin(),g.members.end(),0.0,
            [&](double s,std::size_t id){return s+pockets[id].heuristicScore;})/g.members.size();
        const double interfaceTerm=contacts?edgeSum/contacts:0;
        const double sizePenalty=0.20*std::max(0.0,total/900.0-1.0);
        // A small prior rewards viable site size, not maximally large unions.
        const double sizePrior=0.10*(1-std::exp(-total/300.0));
        g.score=0.55*interfaceTerm+0.35*std::clamp(pocketPrior/0.15,0.0,1.0)
            +sizePrior-sizePenalty;
    };
    auto better=[](const PocketGroup&a,const PocketGroup&b){
        if(a.score!=b.score)return a.score>b.score;
        return a.members<b.members;
    };
    std::vector<PocketGroup> beam,all;
    for(std::size_t i=0;i<n;i++)if(pockets[i].volume<=maxVolume){
        PocketGroup g;g.members={i};score(g);beam.push_back(g);all.push_back(g);
    }
    std::sort(beam.begin(),beam.end(),better);
    if(beam.size()>beamWidth)beam.resize(beamWidth);
    std::set<std::vector<std::size_t>> emitted;
    for(const auto& g:all)emitted.insert(g.members);
    for(std::size_t depth=2;depth<=maxMembers&&!beam.empty();depth++){
        std::vector<PocketGroup> next;
        for(const auto& group:beam)for(auto member:group.members)
            for(const auto& edge:adjacency[member]){
                if(std::binary_search(group.members.begin(),group.members.end(),edge.other))continue;
                PocketGroup g=group;
                g.members.push_back(edge.other);
                std::sort(g.members.begin(),g.members.end());
                if(!emitted.insert(g.members).second)continue;
                score(g);
                if(g.volume>maxVolume||g.score<=-1e8)continue;
                next.push_back(g);all.push_back(g);
            }
        std::sort(next.begin(),next.end(),better);
        if(next.size()>beamWidth)next.resize(beamWidth);
        beam=std::move(next);
    }
    std::sort(all.begin(),all.end(),better);
    if(all.size()>topK)all.resize(topK);
    return all;
}
} // namespace ligint
