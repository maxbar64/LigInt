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
std::vector<AlphaPocket> findAlphaPockets(const std::vector<PocketAtom>& atoms,const AlphaPocketOptions& o){
    if(o.minRadius<=0||o.maxRadius<=o.minRadius||o.neighborCutoff<=0||
       o.clusterDistance<=0||o.minSpheres==0||o.maxCandidates==0)
        throw std::invalid_argument("Invalid alpha pocket options");
    std::vector<AlphaSphere> candidates;
    // Neighbor-limited tetrahedra avoid cubic scanning of the whole structure.
    // The sphere must be empty in the VdW-adjusted sense; atom centers define
    // the circumcenter candidates, with the smallest atom clearance as radius.
    constexpr std::size_t nearestLimit=10;
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
                bool duplicate=false;
                for(const auto& s:candidates)
                    if(dist2(center,s.center)<0.25){duplicate=true;break;}
                if(!duplicate)candidates.push_back({center,radius});
                if(candidates.size()>=o.maxCandidates)break;
            }
        if(candidates.size()>=o.maxCandidates)break;
    }
    DSU groups(candidates.size());
    for(std::size_t i=0;i<candidates.size();++i)
        for(std::size_t j=i+1;j<candidates.size();++j)
            if(dist2(candidates[i].center,candidates[j].center)<=sq(o.clusterDistance))
                groups.join(i,j);
    std::vector<AlphaPocket> result;
    std::vector<std::size_t> roots;
    for(std::size_t i=0;i<candidates.size();++i){
        const auto root=groups.find(i);
        auto it=std::find(roots.begin(),roots.end(),root);
        if(it==roots.end()){roots.push_back(root);result.emplace_back();it=roots.end()-1;}
        auto& pocket=result[static_cast<std::size_t>(it-roots.begin())];
        pocket.spheres.push_back(candidates[i]);
        for(int k=0;k<3;++k)pocket.center[k]+=candidates[i].center[k];
    }
    result.erase(std::remove_if(result.begin(),result.end(),[&](AlphaPocket& p){
        if(p.spheres.size()<o.minSpheres)return true;
        for(int k=0;k<3;++k)p.center[k]/=p.spheres.size();
        for(const auto& s:p.spheres)p.volume+=4.1887902047863905*s.radius*s.radius*s.radius;
        return false;
    }),result.end());
    std::sort(result.begin(),result.end(),[](const AlphaPocket&a,const AlphaPocket&b){
        return a.spheres.size()>b.spheres.size();
    });
    return result;
}
