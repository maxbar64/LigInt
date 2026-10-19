#include "ligint/pocket.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
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
         o.directions>=6 && o.directions<=4096 && o.maxGridPoints>0))
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
    std::sort(out.pockets.begin(),out.pockets.end(),[](const PocketCandidate& a,const PocketCandidate& b){
        if(a.volume!=b.volume)return a.volume>b.volume;
        return a.peakBuriedness>b.peakBuriedness;
    });
    return out;
}
} // namespace ligint
