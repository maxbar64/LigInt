#include "ligint/pocket_alpha.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace ligint;
int main(int argc,char**argv){
 try{
  if(argc<2){std::cerr<<"Usage: ligint-pocket-alpha structure.pdb [--min-radius A] [--max-radius A] [--neighbor-cutoff A] [--cluster-distance A] [--max-cluster-radius A] [--min-contact-atoms N] [--min-spheres N] [--output CSV] [--pdb FILE]\n";return 1;}
  AlphaPocketOptions o;
  std::string input,output,pdb,audit;
  for(int i=1;i<argc;++i){
   const std::string key=argv[i];
   if(key=="--help"||key=="-h"){std::cout<<"Native experimental alpha-sphere-like pocket finder (not fpocket).\nUsage: ligint-pocket-alpha structure.pdb [--min-radius A] [--max-radius A] [--neighbor-cutoff A] [--cluster-distance A] [--min-spheres N] [--output CSV] [--pdb FILE]\n";return 0;}
   if(key.rfind("--",0)!=0){if(!input.empty())throw std::invalid_argument("Only one input structure");input=key;continue;}
   if(i+1>=argc)throw std::invalid_argument("Missing value for "+key);
   const std::string value=argv[++i];
   if(key=="--min-radius")o.minRadius=std::stod(value);
   else if(key=="--max-radius")o.maxRadius=std::stod(value);
   else if(key=="--neighbor-cutoff")o.neighborCutoff=std::stod(value);
   else if(key=="--cluster-distance")o.clusterDistance=std::stod(value);
   else if(key=="--max-cluster-radius")o.maxClusterRadius=std::stod(value);
   else if(key=="--min-contact-atoms")o.minContactAtoms=std::stoull(value);
   else if(key=="--audit")audit=value;
   else if(key=="--min-spheres")o.minSpheres=std::stoull(value);
   else if(key=="--output")output=value;
   else if(key=="--pdb")pdb=value;
   else throw std::invalid_argument("Unknown option "+key);
  }
  if(input.empty())throw std::invalid_argument("Missing PDB input");
  const auto atoms=readPocketAtoms(input);
  const auto pockets=findAlphaPockets(atoms,o);
  std::ofstream file;
  if(!output.empty()){file.open(output);if(!file)throw std::runtime_error("Unable to open CSV output");}
  std::ostream& out=output.empty()?std::cout:file;
  out<<"rank,center_x,center_y,center_z,sphere_count,sphere_volume_sum_A3\n";
  out<<std::fixed<<std::setprecision(3);
  for(std::size_t i=0;i<pockets.size();++i){
   const auto& p=pockets[i];
   out<<i+1<<','<<p.center[0]<<','<<p.center[1]<<','<<p.center[2]<<','
      <<p.spheres.size()<<','<<p.volume<<'\n';
  }
  if(!audit.empty()){
   std::ofstream f(audit);if(!f)throw std::runtime_error("Unable to open audit CSV");
   f<<"rank,sphere_count,min_radius_A,median_radius_A,max_radius_A,max_center_extent_A,volume_sum_A3\n";
   for(std::size_t i=0;i<pockets.size();++i){
    const auto& p=pockets[i];
    std::vector<double> radii;
    double extent=0;
    for(const auto& sp:p.spheres){
     radii.push_back(sp.radius);
     const double dx=sp.center[0]-p.center[0],dy=sp.center[1]-p.center[1],dz=sp.center[2]-p.center[2];
     extent=std::max(extent,std::sqrt(dx*dx+dy*dy+dz*dz));
    }
    std::sort(radii.begin(),radii.end());
    f<<i+1<<','<<radii.size()<<','<<radii.front()<<','<<radii[radii.size()/2]<<','
     <<radii.back()<<','<<extent<<','<<p.volume<<'\n';
   }
  }
  if(!pdb.empty()){
   std::ofstream f(pdb);if(!f)throw std::runtime_error("Unable to open PDB output");
   std::size_t serial=0;
   for(std::size_t i=0;i<pockets.size();++i)for(const auto& sp:pockets[i].spheres){
     if(++serial>99999)break;
     char line[160];
     std::snprintf(line,sizeof(line),"HETATM%5zu  C   SPH A%4zu    %8.3f%8.3f%8.3f  1.00%6.2f           C\n",serial,(i+1)%10000,sp.center[0],sp.center[1],sp.center[2],sp.radius);
     f<<line;
   }
   f<<"END\n";
  }
  return 0;
 }catch(const std::exception&e){std::cerr<<"ligint-pocket-alpha: "<<e.what()<<'\n';return 2;}
}
