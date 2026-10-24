#include "ligint/pocket_alpha.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
namespace ligint {
namespace {
double sq(double x){return x*x;}
double dist2(const V3&a,const V3&b){return sq(a[0]-b[0])+sq(a[1]-b[1])+sq(a[2]-b[2]);}
bool circumcenter(const V3&a,const V3&b,const V3&c,const V3&d,V3&out){
    double m[3][4]{};
    const V3* p[]={&b,&c,&d};
    for(int i=0;i<3;++i){
        double norm=0;
        for(int j=0;j<3;++j){m[i][j]=2*((*p[i])[j]-a[j]);norm+=sq((*p[i])[j])-sq(a[j]);}
        m[i][3]=norm;
    }
    for(int j=0;j<3;++j){
        int pivot=j;
        for(int i=j+1;i<3;++i)if(std::abs(m[i][j])>std::abs(m[pivot][j]))pivot=i;
        if(std::abs(m[pivot][j])<1e-8)return false;
        for(int k=j;k<4;++k)std::swap(m[j][k],m[pivot][k]);
        const double scale=m[j][j];
        for(int k=j;k<4;++k)m[j][k]/=scale;
        for(int i=0;i<3;++i)if(i!=j){
            const double ratio=m[i][j];
            for(int k=j;k<4;++k)m[i][k]-=ratio*m[j][k];
        }
    }
    for(int j=0;j<3;++j)out[j]=m[j][3];
    return true;
}
struct DSU{
 std::vector<std::size_t> parent;
 explicit DSU(std::size_t n):parent(n){std::iota(parent.begin(),parent.end(),0);}
 std::size_t find(std::size_t i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
 void join(std::size_t i,std::size_t j){i=find(i);j=find(j);if(i!=j)parent[j]=i;}
};
}
std::vector<AlphaPocket> findAlphaPockets(const std::vector<PocketAtom>& atoms,const AlphaPocketOptions& o, AlphaMergeAudit* audit){
    if(o.minRadius<=0||o.maxRadius<=o.minRadius||o.neighborCutoff<=0||
       o.clusterDistance<=0||o.maxClusterRadius<=0||o.macroLinkDistance<=0||o.macroMaxDiameter<=0||o.macroMinClearance<0||o.minContactAtoms<4||o.minSpheres==0||o.maxCandidates==0)
        throw std::invalid_argument("Invalid alpha pocket options");
    std::vector<AlphaSphere> candidates;
    // Neighbor-limited tetrahedra avoid cubic scanning of the whole structure.
    // The sphere must be empty in the VdW-adjusted sense; atom centers define
    // the circumcenter candidates, with the smallest atom clearance as radius.
    constexpr std::size_t nearestLimit=7;
    for(std::size_t i=0;i<atoms.size();++i){
        std::vector<std::pair<double,std::size_t>> neighbors;
        for(std::size_t j=0;j<atoms.size();++j)if(i!=j){
            const double d=dist2(atoms[i].pos,atoms[j].pos);
            if(d<=sq(o.neighborCutoff))neighbors.push_back({d,j});
        }
        if(neighbors.size()<3)continue;
        const std::size_t n=std::min(nearestLimit,neighbors.size());
        std::partial_sort(neighbors.begin(),neighbors.begin()+n,neighbors.end());
        for(std::size_t a=0;a<n;++a)for(std::size_t b=a+1;b<n;++b)
            for(std::size_t c=b+1;c<n;++c){
                V3 center{};
                if(!circumcenter(atoms[i].pos,atoms[neighbors[a].second].pos,
                    atoms[neighbors[b].second].pos,atoms[neighbors[c].second].pos,center))continue;
                // Local tetrahedral atoms identify a candidate Voronoi position.
                // Use the minimum clearance from all atom VdW surfaces as the
                // actual sphere radius, including mixed atomic radii.
                double radius=o.maxRadius+1;
                for(const auto& atom:atoms)
                    radius=std::min(radius,std::sqrt(dist2(center,atom.pos))-atom.radius);
                if(radius<o.minRadius||radius>o.maxRadius)continue;
                std::size_t contacts=0;
                for(const auto& atom:atoms)
                    if(std::sqrt(dist2(center,atom.pos))-atom.radius<=radius+2.0)++contacts;
                if(contacts<o.minContactAtoms)continue;
                bool duplicate=false;
                for(const auto& s:candidates)
                    if(dist2(center,s.center)<0.25){duplicate=true;break;}
                if(!duplicate)candidates.push_back({center,radius});
                if(candidates.size()>=o.maxCandidates)break;
            }
        if(candidates.size()>=o.maxCandidates)break;
    }
    // Bounded centroid clustering: no unrestricted single-link chaining.
    std::sort(candidates.begin(),candidates.end(),[](const auto&a,const auto&b){return a.radius>b.radius;});
    std::vector<AlphaPocket> result;
    for(const auto& sphere:candidates){
        std::size_t best=result.size();
        double bestDistance=sq(o.clusterDistance);
        for(std::size_t i=0;i<result.size();++i){
            const double d2=dist2(sphere.center,result[i].center);
            if(d2>=bestDistance)continue;
            bool bounded=true;
            for(const auto& member:result[i].spheres)
                if(dist2(sphere.center,member.center)>sq(o.maxClusterRadius*2)){bounded=false;break;}
            if(bounded){best=i;bestDistance=d2;}
        }
        if(best==result.size()){
            AlphaPocket p;p.center=sphere.center;p.spheres.push_back(sphere);
            result.push_back(std::move(p));
        }else{
            auto& p=result[best];
            const double n=static_cast<double>(p.spheres.size());
            for(int k=0;k<3;++k)p.center[k]=(p.center[k]*n+sphere.center[k])/(n+1);
            p.spheres.push_back(sphere);
        }
    }
    result.erase(std::remove_if(result.begin(),result.end(),[&](AlphaPocket& p){
        if(p.spheres.size()<o.minSpheres)return true;
        for(const auto& s:p.spheres)p.volume+=4.1887902047863905*s.radius*s.radius*s.radius;
        return false;
    }),result.end());
    if(audit){audit->microclusters=result.size();audit->decisions.clear();}
    if(o.mergeMicroclusters && result.size()>1) {
        // Evaluate candidate microcluster links from the closest pair of sphere
        // centers, not solely centroid distances. Merge only if the path is free
        // of atom VdW surfaces and the resulting center-span remains bounded.
        struct Edge {std::size_t a,b;double distance;V3 from,to;};
        std::vector<Edge> edges;
        for(std::size_t i=0;i<result.size();++i)for(std::size_t j=i+1;j<result.size();++j){
            double best=sq(o.macroLinkDistance);
            V3 first{},second{};
            bool found=false;
            for(const auto& x:result[i].spheres)for(const auto& y:result[j].spheres){
                const double d=dist2(x.center,y.center);
                if(d<best){best=d;first=x.center;second=y.center;found=true;}
            }
            if(found)edges.push_back({i,j,std::sqrt(best),first,second});
        }
        std::sort(edges.begin(),edges.end(),[](const Edge& a,const Edge& b){return a.distance<b.distance;});
        DSU components(result.size());
        for(const auto& e:edges){
            const auto a=components.find(e.a),b=components.find(e.b);
            if(a==b)continue;
            AlphaMergeDecision decision;
            decision.first=e.a+1;decision.second=e.b+1;decision.linkDistance=e.distance;
            // Geometric accessibility along the candidate bridge, sampled
            // at 0.5 Å intervals. This is not a full channel traversal.
            double minClearance=1e30;
            const int steps=std::max(1,static_cast<int>(std::ceil(e.distance/0.5)));
            for(int step=0;step<=steps;++step){
                const double t=static_cast<double>(step)/steps;
                V3 p{};
                for(int k=0;k<3;++k)p[k]=e.from[k]*(1-t)+e.to[k]*t;
                for(const auto& atom:atoms)
                    minClearance=std::min(minClearance,std::sqrt(dist2(p,atom.pos))-atom.radius);
            }
            decision.minClearance=minClearance;
            decision.endpointClearance=1e30;
            for(const auto& atom:atoms) {
                decision.endpointClearance=std::min(decision.endpointClearance,
                    std::sqrt(dist2(e.from,atom.pos))-atom.radius);
                decision.endpointClearance=std::min(decision.endpointClearance,
                    std::sqrt(dist2(e.to,atom.pos))-atom.radius);
            }
            if(minClearance<o.macroMinClearance){
                decision.reason="insufficient_clearance";
            }else{
                double diameter=0;
                // Bound the span of ALL microcluster centers in the proposed
                // union: prevents transitive chaining through many short links.
                for(std::size_t i=0;i<result.size();++i){
                    const auto ri=components.find(i);
                    if(ri!=a && ri!=b)continue;
                    for(std::size_t j=i+1;j<result.size();++j){
                        const auto rj=components.find(j);
                        if(rj!=a && rj!=b)continue;
                        diameter=std::max(diameter,std::sqrt(dist2(result[i].center,result[j].center)));
                    }
                }
                decision.mergedDiameter=diameter;
                if(diameter>o.macroMaxDiameter)decision.reason="diameter_limit";
                else {decision.accepted=true;decision.reason="accepted";components.join(a,b);}
            }
            if(audit)audit->decisions.push_back(decision);
        }
        std::vector<AlphaPocket> macro;
        std::vector<std::size_t> labels;
        for(std::size_t i=0;i<result.size();++i){
            const auto root=components.find(i);
            auto it=std::find(labels.begin(),labels.end(),root);
            if(it==labels.end()){labels.push_back(root);macro.emplace_back();it=labels.end()-1;}
            auto& dst=macro[static_cast<std::size_t>(it-labels.begin())];
            const auto& src=result[i];
            for(const auto& sphere:src.spheres)dst.spheres.push_back(sphere);
            dst.volume+=src.volume;
            for(int k=0;k<3;++k)dst.center[k]+=src.center[k]*src.spheres.size();
        }
        for(auto& pocket:macro)
            for(int k=0;k<3;++k)pocket.center[k]/=pocket.spheres.size();
        result=std::move(macro);
    }
    std::sort(result.begin(),result.end(),[&](const AlphaPocket&a,const AlphaPocket&b){
        if(o.ranking==AlphaPocketRanking::SphereVolume && a.volume!=b.volume)
            return a.volume>b.volume;
        if(a.spheres.size()!=b.spheres.size())
            return a.spheres.size()>b.spheres.size();
        return a.volume>b.volume;
    });
    return result;
}

} // namespace ligint
