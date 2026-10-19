#pragma once
#include "ligint/pocket.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
namespace gui {
struct SmoothTriangle {std::array<ligint::V3,3> vertices;};
inline std::vector<SmoothTriangle> makeSmoothPocketMesh(
    const std::vector<ligint::PocketPoint>& points,double spacing,int passes=2){
    std::vector<SmoothTriangle> triangles;
    if(points.empty()||spacing<=0)return triangles;
    ligint::V3 low=points[0].pos,high=low;
    for(const auto& p:points)for(int k=0;k<3;k++){
        low[k]=std::min(low[k],p.pos[k]);high[k]=std::max(high[k],p.pos[k]);
    }
    std::array<int,3> n;
    for(int k=0;k<3;k++){low[k]-=2*spacing;n[k]=static_cast<int>(std::ceil((high[k]-low[k])/spacing))+5;}
    const std::size_t count=static_cast<std::size_t>(n[0])*n[1]*n[2];
    if(count>3500000)return triangles; // safeguard large selections
    auto idx=[&](int x,int y,int z){return (static_cast<std::size_t>(z)*n[1]+y)*n[0]+x;};
    std::vector<float> field(count,0),next(count,0);
    for(const auto& p:points){
        int xyz[3];
        for(int k=0;k<3;k++)xyz[k]=static_cast<int>(std::llround((p.pos[k]-low[k])/spacing));
        if(xyz[0]>=0&&xyz[0]<n[0]&&xyz[1]>=0&&xyz[1]<n[1]&&xyz[2]>=0&&xyz[2]<n[2])
            field[idx(xyz[0],xyz[1],xyz[2])]=1;
    }
    passes=std::clamp(passes,0,4);
    for(int pass=0;pass<passes;pass++){
        for(int z=1;z<n[2]-1;z++)for(int y=1;y<n[1]-1;y++)for(int x=1;x<n[0]-1;x++){
            float sum=0;
            for(int dz=-1;dz<=1;dz++)for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){
                const int w=(dx==0?2:1)*(dy==0?2:1)*(dz==0?2:1);
                sum+=w*field[idx(x+dx,y+dy,z+dz)];
            }
            next[idx(x,y,z)]=sum/64.0f;
        }
        field.swap(next);
        std::fill(next.begin(),next.end(),0);
    }
    const float iso=passes?0.35f:0.5f;
    constexpr int corners[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    constexpr int tetra[6][4]={{0,5,1,6},{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6}};
    auto interp=[&](const ligint::V3& a,const ligint::V3& b,float fa,float fb){
        const double t=std::clamp(double(iso-fa)/double(fb-fa),0.0,1.0);
        ligint::V3 v;for(int k=0;k<3;k++)v[k]=a[k]+t*(b[k]-a[k]);return v;
    };
    for(int z=0;z<n[2]-1;z++)for(int y=0;y<n[1]-1;y++)for(int x=0;x<n[0]-1;x++){
        ligint::V3 v[8];float f[8];bool anyIn=false,anyOut=false;
        for(int c=0;c<8;c++){
            const int ix=x+corners[c][0],iy=y+corners[c][1],iz=z+corners[c][2];
            v[c]={low[0]+ix*spacing,low[1]+iy*spacing,low[2]+iz*spacing};
            f[c]=field[idx(ix,iy,iz)];anyIn|=f[c]>=iso;anyOut|=f[c]<iso;
        }
        if(!anyIn||!anyOut)continue;
        for(const auto& t:tetra){
            int in[4],out[4],ni=0,no=0;
            for(int j=0;j<4;j++)(f[t[j]]>=iso?in[ni++]:out[no++])=t[j];
            if(ni==0||ni==4)continue;
            if(ni==1||ni==3){
                const int single=ni==1?in[0]:out[0];
                const int* other=ni==1?out:in;
                triangles.push_back({{{interp(v[single],v[other[0]],f[single],f[other[0]]),
                                       interp(v[single],v[other[1]],f[single],f[other[1]]),
                                       interp(v[single],v[other[2]],f[single],f[other[2]])}}});
            }else{
                const auto a=interp(v[in[0]],v[out[0]],f[in[0]],f[out[0]]);
                const auto b=interp(v[in[0]],v[out[1]],f[in[0]],f[out[1]]);
                const auto c=interp(v[in[1]],v[out[0]],f[in[1]],f[out[0]]);
                const auto d=interp(v[in[1]],v[out[1]],f[in[1]],f[out[1]]);
                triangles.push_back({{{a,b,c}}});triangles.push_back({{{b,d,c}}});
            }
        }
    }
    return triangles;
}
}
