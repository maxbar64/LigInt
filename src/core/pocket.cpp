#include "ligint/pocket.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <queue>
#include <stdexcept>

namespace ligint {
namespace {
struct Grid {
    V3 origin{0,0,0};
    int nx=0,ny=0,nz=0;
    double step=1;
    std::size_t index(int x,int y,int z) const {
        return (static_cast<std::size_t>(z)*ny+y)*nx+x;
    }
    bool inside(int x,int y,int z) const {
        return x>=0 && y>=0 && z>=0 && x<nx && y<ny && z<nz;
    }
    V3 coord(int x,int y,int z) const {
        return {origin[0]+step*x,origin[1]+step*y,origin[2]+step*z};
    }
};
std::vector<V3> sphereDirections(int n) {
    std::vector<V3> dirs;
    dirs.reserve(static_cast<std::size_t>(n));
    constexpr double pi=3.14159265358979323846;
    const double golden=pi*(3.0-std::sqrt(5.0));
    for(int k=0;k<n;k++) {
        double y=1.0-2.0*(k+0.5)/n;
        double r=std::sqrt(std::max(0.0,1.0-y*y));
        double t=k*golden;
        dirs.push_back({r*std::cos(t),y,r*std::sin(t)});
    }
    return dirs;
}
}
PocketSearchResult findPockets(const std::vector<PocketAtom>& input, const PocketOptions& o) {
    if(!(std::isfinite(o.spacing) && o.spacing>0 && std::isfinite(o.rayLength) && o.rayLength>0 &&
         std::isfinite(o.maxSurfaceDistance) && o.maxSurfaceDistance>0 &&
         std::isfinite(o.minBuriedness) && o.minBuriedness>=0 && o.minBuriedness<=1 &&
         o.directions>=6 && o.directions<=4096 && o.maxGridPoints>0 &&
         std::isfinite(o.seedSpacing) && o.seedSpacing>0 &&
         std::isfinite(o.seedProminence) && o.seedProminence>=0 && o.seedProminence<=1 && o.maxSeeds>0 &&
         std::isfinite(o.mergeSaddleDrop) && o.mergeSaddleDrop>=0 && o.mergeSaddleDrop<=1 &&
         std::isfinite(o.mergeMaxSeedDistance) && o.mergeMaxSeedDistance>0 &&
         std::isfinite(o.mergeMinContactArea) && o.mergeMinContactArea>=0))
        throw std::invalid_argument("Invalid pocket search options");
    std::vector<PocketAtom> atoms;
    for(const auto& a:input) {
        if(!o.chains.empty() && std::find(o.chains.begin(),o.chains.end(),a.chain)==o.chains.end()) continue;
        if(!std::isfinite(a.radius) || a.radius<=0) continue;
        bool valid=true;
        for(double v:a.pos) if(!std::isfinite(v)) valid=false;
        if(valid) atoms.push_back(a);
    }
    PocketSearchResult out;
    out.selectedAtoms=atoms.size();
    if(atoms.empty()) return out;
    V3 lo=atoms.front().pos,hi=lo;
    double maxRadius=0;
    for(const auto& a:atoms) {
        maxRadius=std::max(maxRadius,a.radius);
        for(int d=0;d<3;d++){lo[d]=std::min(lo[d],a.pos[d]);hi[d]=std::max(hi[d],a.pos[d]);}
    }
    Grid g;
    g.step=o.spacing;
    const double margin=maxRadius+o.maxSurfaceDistance+o.spacing;
    int dims[3];
    for(int d=0;d<3;d++) {
        g.origin[d]=lo[d]-margin;
        double len=hi[d]-lo[d]+2*margin;
        if(len/o.spacing>1000000) throw std::length_error("Pocket grid too large");
        dims[d]=static_cast<int>(std::ceil(len/o.spacing))+1;
    }
    g.nx=dims[0];g.ny=dims[1];g.nz=dims[2];
    const std::size_t xy=static_cast<std::size_t>(g.nx)*g.ny;
    if(xy>o.maxGridPoints/static_cast<std::size_t>(g.nz))
        throw std::length_error("Pocket grid exceeds maxGridPoints");
    const std::size_t count=xy*g.nz;
    // Occupied grid sites approximate a van der Waals union.
    std::vector<std::uint8_t> occupied(count,0);
    for(const auto& a:atoms) {
        int low[3],high[3];
        for(int d=0;d<3;d++) {
            low[d]=std::max(0,static_cast<int>(std::floor((a.pos[d]-a.radius-g.origin[d])/g.step)));
            high[d]=std::min(dims[d]-1,static_cast<int>(std::ceil((a.pos[d]+a.radius-g.origin[d])/g.step)));
        }
        for(int z=low[2];z<=high[2];z++)for(int y=low[1];y<=high[1];y++)for(int x=low[0];x<=high[0];x++){
            const V3 p=g.coord(x,y,z);
            const double dx=p[0]-a.pos[0],dy=p[1]-a.pos[1],dz=p[2]-a.pos[2];
            if(dx*dx+dy*dy+dz*dz<=a.radius*a.radius)occupied[g.index(x,y,z)]=1;
        }
    }
    // Multi-source distance in 26-neighborhood: candidate pruning only.
    // Dijkstra-like distances are not required for buriedness: keep a conservative
    // Euclidean cutoff using atom centers for each sampled voxel instead.
    std::vector<float> score(count,-1);
    const auto dirs=sphereDirections(o.directions);
    const double rayStep=std::min(o.spacing*0.5,0.5);
    const int nsteps=static_cast<int>(std::ceil(o.rayLength/rayStep));
    for(int z=0;z<g.nz;z++)for(int y=0;y<g.ny;y++)for(int x=0;x<g.nx;x++){
        const std::size_t id=g.index(x,y,z);
        if(occupied[id])continue;
        const V3 p=g.coord(x,y,z);
        bool near=false;
        for(const auto& a:atoms){
            double dx=p[0]-a.pos[0],dy=p[1]-a.pos[1],dz=p[2]-a.pos[2];
            double lim=a.radius+o.maxSurfaceDistance;
            if(dx*dx+dy*dy+dz*dz<=lim*lim){near=true;break;}
        }
        if(!near)continue;
        out.sampledPoints++;
        int blocked=0;
        for(const auto& u:dirs){
            bool hit=false;
            for(int t=1;t<=nsteps;t++){
                double s=std::min(t*rayStep,o.rayLength);
                int xx=static_cast<int>(std::llround((p[0]+s*u[0]-g.origin[0])/g.step));
                int yy=static_cast<int>(std::llround((p[1]+s*u[1]-g.origin[1])/g.step));
                int zz=static_cast<int>(std::llround((p[2]+s*u[2]-g.origin[2])/g.step));
                if(!g.inside(xx,yy,zz))break;
                if(occupied[g.index(xx,yy,zz)]){hit=true;break;}
            }
            if(hit)blocked++;
        }
        float b=static_cast<float>(blocked)/static_cast<float>(dirs.size());
        if(b>=o.minBuriedness)score[id]=b;
    }
    std::vector<std::uint8_t> seen(count,0);
    // 6-connected components; subsequent refinement can split wide pockets.
    for(int z=0;z<g.nz;z++)for(int y=0;y<g.ny;y++)for(int x=0;x<g.nx;x++){
        const auto start=g.index(x,y,z);
        if(seen[start] || score[start]<0)continue;
        seen[start]=1;
        std::deque<std::array<int,3>> queue;
        queue.push_back({x,y,z});
        PocketCandidate pocket;
        while(!queue.empty()){
            auto q=queue.front();queue.pop_front();
            const auto id=g.index(q[0],q[1],q[2]);
            const V3 p=g.coord(q[0],q[1],q[2]);
            pocket.points.push_back({p,score[id]});
            pocket.peakBuriedness=std::max(pocket.peakBuriedness,static_cast<double>(score[id]));
            for(int d=0;d<3;d++)pocket.center[d]+=p[d];
            const int dx[6]={1,-1,0,0,0,0};
            const int dy[6]={0,0,1,-1,0,0};
            const int dz[6]={0,0,0,0,1,-1};
            for(int k=0;k<6;k++){
                int xx=q[0]+dx[k],yy=q[1]+dy[k],zz=q[2]+dz[k];
                if(!g.inside(xx,yy,zz))continue;
                auto j=g.index(xx,yy,zz);
                if(!seen[j] && score[j]>=0){seen[j]=1;queue.push_back({xx,yy,zz});}
            }
        }
        pocket.pointCount=pocket.points.size();
        if(pocket.pointCount<o.minClusterPoints)continue;
        for(auto& c:pocket.center)c/=static_cast<double>(pocket.pointCount);
        pocket.volume=pocket.pointCount*g.step*g.step*g.step;
        out.pockets.push_back(std::move(pocket));
    }
    // Seeded 3D watershed on the precomputed buriedness grid. Smooth only for seed
    // selection/propagation: keep the original buriedness on exported pseudoatoms.
    if(o.segment) {
        std::vector<PocketCandidate> separated;
        for(const auto& connected:out.pockets) {
            if(connected.points.size()<2*o.minClusterPoints) {
                separated.push_back(connected);
                continue;
            }
            std::vector<std::size_t> ids;
            ids.reserve(connected.points.size());
            for(const auto& pt:connected.points) {
                int x=static_cast<int>(std::llround((pt.pos[0]-g.origin[0])/g.step));
                int y=static_cast<int>(std::llround((pt.pos[1]-g.origin[1])/g.step));
                int z=static_cast<int>(std::llround((pt.pos[2]-g.origin[2])/g.step));
                ids.push_back(g.index(x,y,z));
            }
            // Six-neighbor averaging is deterministic and removes single-voxel peaks.
            std::vector<float> smooth(ids.size(),0);
            const int dx[6]={1,-1,0,0,0,0},dy[6]={0,0,1,-1,0,0},dz[6]={0,0,0,0,1,-1};
            for(std::size_t i=0;i<ids.size();i++) {
                auto id=ids[i];
                int x=static_cast<int>(id%g.nx);
                int y=static_cast<int>((id/g.nx)%g.ny);
                int z=static_cast<int>(id/(static_cast<std::size_t>(g.nx)*g.ny));
                double sum=score[id]*2.0,weight=2.0;
                for(int k=0;k<6;k++){
                    if(!g.inside(x+dx[k],y+dy[k],z+dz[k]))continue;
                    const auto j=g.index(x+dx[k],y+dy[k],z+dz[k]);
                    if(score[j]<0)continue;
                    sum+=score[j];weight+=1;
                }
                smooth[i]=static_cast<float>(sum/weight);
            }
            std::vector<std::size_t> order(ids.size());
            for(std::size_t i=0;i<order.size();i++)order[i]=i;
            std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){
                if(smooth[a]!=smooth[b])return smooth[a]>smooth[b];
                return ids[a]<ids[b];
            });
            std::vector<std::size_t> seeds;
            const double minDist2=o.seedSpacing*o.seedSpacing;
            // Greedy nonmaximum suppression. A separation radius avoids a seed
            // for each voxel in a broad plateau with buriedness exactly 1.
            for(auto i:order){
                if(smooth[i]<o.minBuriedness+o.seedProminence)break;
                bool near=false;
                for(auto s:seeds){
                    const auto& a=connected.points[i].pos;
                    const auto& b=connected.points[s].pos;
                    double d2=0;
                    for(int k=0;k<3;k++){double d=a[k]-b[k];d2+=d*d;}
                    if(d2<minDist2){near=true;break;}
                }
                if(!near)seeds.push_back(i);
                if(seeds.size()>=o.maxSeeds)break;
            }
            if(seeds.size()<2){separated.push_back(connected);continue;}
            std::vector<int> belongs(count,-1);
            for(std::size_t i=0;i<ids.size();i++)belongs[ids[i]]=static_cast<int>(i);
            std::vector<int> labels(ids.size(),-1);
            struct Node {float score; std::size_t idx; int label;};
            auto cmp=[](const Node& a,const Node& b){
                if(a.score!=b.score)return a.score<b.score;
                if(a.label!=b.label)return a.label>b.label;
                return a.idx>b.idx;
            };
            std::priority_queue<Node,std::vector<Node>,decltype(cmp)> frontier(cmp);
            for(std::size_t k=0;k<seeds.size();k++)
                frontier.push({smooth[seeds[k]],seeds[k],static_cast<int>(k)});
            while(!frontier.empty()){
                auto node=frontier.top();frontier.pop();
                if(labels[node.idx]!=-1)continue;
                labels[node.idx]=node.label;
                auto id=ids[node.idx];
                int x=static_cast<int>(id%g.nx);
                int y=static_cast<int>((id/g.nx)%g.ny);
                int z=static_cast<int>(id/(static_cast<std::size_t>(g.nx)*g.ny));
                for(int k=0;k<6;k++){
                    if(!g.inside(x+dx[k],y+dy[k],z+dz[k]))continue;
                    int j=belongs[g.index(x+dx[k],y+dy[k],z+dz[k])];
                    if(j>=0 && labels[static_cast<std::size_t>(j)]<0)
                        frontier.push({std::min(node.score,smooth[static_cast<std::size_t>(j)]),
                                       static_cast<std::size_t>(j),node.label});
                }
            }
            // Measure the strongest shared saddle between adjacent watershed basins.
            // Merge only when the saddle is near the weaker basin's maximum.
            if(o.mergeBasins) {
                const std::size_t ns=seeds.size();
                std::vector<float> peak(ns,0.f);
                std::vector<std::size_t> sizes(ns,0);
                for(std::size_t i=0;i<labels.size();i++)if(labels[i]>=0){
                    const auto k=static_cast<std::size_t>(labels[i]);
                    peak[k]=std::max(peak[k],smooth[i]);
                    ++sizes[k];
                }
                std::vector<float> saddle(ns*ns,-1.f);
                std::vector<std::size_t> contact(ns*ns,0);
                for(std::size_t i=0;i<ids.size();i++){
                    const int a=labels[i];
                    if(a<0)continue;
                    const auto id=ids[i];
                    const int x=static_cast<int>(id%g.nx);
                    const int y=static_cast<int>((id/g.nx)%g.ny);
                    const int z=static_cast<int>(id/(static_cast<std::size_t>(g.nx)*g.ny));
                    for(int k=0;k<6;k+=2){
                        const int xx=x+dx[k],yy=y+dy[k],zz=z+dz[k];
                        if(!g.inside(xx,yy,zz))continue;
                        const int j=belongs[g.index(xx,yy,zz)];
                        if(j<0)continue;
                        const int b=labels[static_cast<std::size_t>(j)];
                        if(b<0||a==b)continue;
                        const float v=std::min(smooth[i],smooth[static_cast<std::size_t>(j)]);
                        ++contact[static_cast<std::size_t>(a)*ns+static_cast<std::size_t>(b)];
                        ++contact[static_cast<std::size_t>(b)*ns+static_cast<std::size_t>(a)];
                        saddle[static_cast<std::size_t>(a)*ns+static_cast<std::size_t>(b)]=
                            std::max(saddle[static_cast<std::size_t>(a)*ns+static_cast<std::size_t>(b)],v);
                        saddle[static_cast<std::size_t>(b)*ns+static_cast<std::size_t>(a)]=
                            std::max(saddle[static_cast<std::size_t>(b)*ns+static_cast<std::size_t>(a)],v);
                    }
                }
                std::vector<std::size_t> parent(ns);
                std::vector<std::vector<std::size_t>> members(ns);
                for(std::size_t i=0;i<ns;i++)members[i].push_back(i);
                for(std::size_t i=0;i<ns;i++)parent[i]=i;
                auto root=[&](std::size_t x){while(parent[x]!=x)x=parent[x]=parent[parent[x]];return x;};
                struct Edge{std::size_t a,b;float level;};
                std::vector<Edge> edges;
                for(std::size_t a=0;a<ns;a++)for(std::size_t b=a+1;b<ns;b++){
                    const float s=saddle[a*ns+b];
                    if(s>=0)edges.push_back({a,b,s});
                }
                std::sort(edges.begin(),edges.end(),[](const Edge& a,const Edge& b){
                    if(a.level!=b.level)return a.level>b.level;
                    if(a.a!=b.a)return a.a<b.a;
                    return a.b<b.b;
                });
                for(const auto& edge:edges) {
                    auto a=root(edge.a),b=root(edge.b);
                    if(a==b)continue;
                    const float weaker=std::min(peak[a],peak[b]);
                    if(weaker-edge.level>o.mergeSaddleDrop)continue;
                    if(static_cast<double>(contact[edge.a*ns+edge.b])*g.step*g.step < o.mergeMinContactArea)continue;
                    bool close=true;
                    const double maxD2=o.mergeMaxSeedDistance*o.mergeMaxSeedDistance;
                    for(auto ia:members[a])for(auto ib:members[b]){
                        const auto& pa=connected.points[seeds[ia]].pos;
                        const auto& pb=connected.points[seeds[ib]].pos;
                        double d2=0;
                        for(int k=0;k<3;k++){const double d=pa[k]-pb[k];d2+=d*d;}
                        if(d2>maxD2)close=false;
                    }
                    if(!close)continue;
                    // The stronger maximum remains the representative seed.
                    if(peak[a]<peak[b] || (peak[a]==peak[b]&&a>b))std::swap(a,b);
                    parent[b]=a;
                    members[a].insert(members[a].end(),members[b].begin(),members[b].end());
                    members[b].clear();
                    sizes[a]+=sizes[b];
                    peak[a]=std::max(peak[a],peak[b]);
                }
                for(auto& label:labels)if(label>=0)
                    label=static_cast<int>(root(static_cast<std::size_t>(label)));
            }
            std::vector<PocketCandidate> basins(seeds.size());
            for(std::size_t i=0;i<ids.size();i++){
                const int label=labels[i];
                if(label<0)continue;
                auto& pocket=basins[static_cast<std::size_t>(label)];
                const auto& pt=connected.points[i];
                pocket.points.push_back(pt);
                for(int k=0;k<3;k++)pocket.center[k]+=pt.pos[k];
                pocket.peakBuriedness=std::max(pocket.peakBuriedness,static_cast<double>(pt.buriedness));
            }
            for(auto& pocket:basins){
                pocket.pointCount=pocket.points.size();
                if(pocket.pointCount<o.minClusterPoints)continue;
                for(auto& v:pocket.center)v/=static_cast<double>(pocket.pointCount);
                pocket.volume=pocket.pointCount*g.step*g.step*g.step;
                separated.push_back(std::move(pocket));
            }
        }
        out.pockets=std::move(separated);
    }
    std::sort(out.pockets.begin(),out.pockets.end(),[](const PocketCandidate& a,const PocketCandidate& b){
        if(a.volume!=b.volume)return a.volume>b.volume;
        return a.peakBuriedness>b.peakBuriedness;
    });
    return out;
}
} // namespace ligint
