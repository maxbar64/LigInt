#include "ligint/pocket.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace ligint;

static void help(){
    std::cout<<"Usage: ligint-pocket structure.pdb|structure.cif [options]\n"
        <<"  --chains A,B        selected protein chains (default all)\n"
        <<"  --spacing A         voxel spacing, default 1.0\n"
        <<"  --ray-length A      ray length, default 10.0\n"
        <<"  --directions N      ray directions, default 42\n"
        <<"  --min-buriedness B  threshold 0..1, default 0.5\n"
        <<"  --surface-distance A maximum candidate distance to atoms, default 5.0\n"
        <<"  --min-points N      minimum cluster size, default 30\n"
        <<"  --top N             output at most N pockets (default all)\n"
        <<"  --rank-by volume|heuristic  pocket ordering (default volume)\n"
        <<"  --no-segmentation   retain original connected components\n"
        <<"  --seed-spacing A    minimum seed separation (default 6 A)\n"
        <<"  --seed-prominence B minimum seed buriedness over cutoff (default 0.10)\n"
        <<"  --max-seeds N       maximum seeds per connected region (default 64)\n"
        <<"  --merge-saddle-drop B maximum buriedness peak-to-saddle drop (default 0.08)\n"
        <<"  --merge-max-distance A maximum seed-group span (default 10 A)\n"
        <<"  --merge-min-contact A2 minimum shared voxel-face area (default 8 A2)\n"
        <<"  --no-basin-merge   disable merging of adjacent watershed basins\n"
        <<"  --output FILE       write summary as JSON (default stdout)\n"
        <<"  --points FILE       export selected buried voxels as CSV\n"
        <<"  --pdb FILE          export pocket pseudoatoms (B = 100 * buriedness)\n"
        <<"  --format auto|pdb|cif\n";
}
static double number(const std::string& value){
    std::size_t used=0;
    double x=std::stod(value,&used);
    if(used!=value.size()||!std::isfinite(x))throw std::invalid_argument("Invalid numeric value: "+value);
    return x;
}
static std::size_t nonnegativeInteger(const std::string& value){
    if(value.empty()||value.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("Invalid integer: "+value);
    return static_cast<std::size_t>(std::stoull(value));
}
int main(int argc,char** argv){
    try{
        if(argc<2){help();return 1;}
        PocketOptions o;
        std::string file,output,points,pdb,rankBy="volume";
        StructureFormat format=StructureFormat::Auto;
        std::size_t top=0;
        for(int i=1;i<argc;i++){
            std::string a=argv[i];
            if(a=="--help" || a=="-h"){help();return 0;}
            if(a.substr(0,2)!="--"){
                if(!file.empty())throw std::invalid_argument("Only one structure file is supported");
                file=a;continue;
            }
            if(a=="--no-segmentation"){o.segment=false;continue;}
            if(a=="--no-basin-merge"){o.mergeBasins=false;continue;}
            if(i+1>=argc)throw std::invalid_argument("Missing value for "+a);
            std::string v=argv[++i];
            if(a=="--chains"){
                std::istringstream s(v);std::string chain;
                while(std::getline(s,chain,',')){
                    if(chain.empty())throw std::invalid_argument("Empty chain ID");
                    o.chains.push_back(chain);
                }
                if(o.chains.empty())throw std::invalid_argument("Empty chain selection");
            }else if(a=="--spacing")o.spacing=number(v);
            else if(a=="--ray-length")o.rayLength=number(v);
            else if(a=="--directions"){
                auto n=nonnegativeInteger(v);
                if(n>4096)throw std::invalid_argument("Too many ray directions");
                o.directions=static_cast<int>(n);
            }else if(a=="--min-buriedness")o.minBuriedness=number(v);
            else if(a=="--surface-distance")o.maxSurfaceDistance=number(v);
            else if(a=="--min-points")o.minClusterPoints=nonnegativeInteger(v);
            else if(a=="--top")top=nonnegativeInteger(v);
            else if(a=="--rank-by"){
                if(v!="volume" && v!="heuristic")throw std::invalid_argument("Invalid ranking: "+v);
                rankBy=v;
            }
            else if(a=="--seed-spacing")o.seedSpacing=number(v);
            else if(a=="--seed-prominence")o.seedProminence=number(v);
            else if(a=="--max-seeds")o.maxSeeds=nonnegativeInteger(v);
            else if(a=="--merge-saddle-drop")o.mergeSaddleDrop=number(v);
            else if(a=="--merge-max-distance")o.mergeMaxSeedDistance=number(v);
            else if(a=="--merge-min-contact")o.mergeMinContactArea=number(v);
            else if(a=="--output")output=v;
            else if(a=="--points")points=v;
            else if(a=="--pdb")pdb=v;
            else if(a=="--format"){
                if(v=="pdb")format=StructureFormat::Pdb;
                else if(v=="cif"||v=="mmcif")format=StructureFormat::Mmcif;
                else if(v=="auto")format=StructureFormat::Auto;
                else throw std::invalid_argument("Invalid structure format");
            }else throw std::invalid_argument("Unknown option: "+a);
        }
        if(file.empty())throw std::invalid_argument("Missing structure file");
        auto atoms=readPocketAtoms(file,format);
        auto r=findPockets(atoms,o);
        if(!o.chains.empty()&&r.selectedAtoms==0)throw std::runtime_error("Selected chains contain no protein atoms");
        if(rankBy=="heuristic")std::stable_sort(r.pockets.begin(),r.pockets.end(),
            [](const PocketCandidate& a,const PocketCandidate& b){return a.heuristicScore>b.heuristicScore;});
        std::ofstream dest,pointDest,pdbDest;
        std::ostream* stream=&std::cout;
        if(!output.empty()){
            dest.open(output);
            if(!dest)throw std::runtime_error("Cannot write: "+output);
            stream=&dest;
        }
        const std::size_t n=top?std::min(top,r.pockets.size()):r.pockets.size();
        *stream<<std::setprecision(7)<<"{\n  \"selected_atoms\": "<<r.selectedAtoms
               <<",\n  \"sampled_points\": "<<r.sampledPoints
               <<",\n  \"pocket_count\": "<<r.pockets.size()<<",\n  \"pockets\": [\n";
        for(std::size_t i=0;i<n;i++){
            const auto& p=r.pockets[i];
            *stream<<"    {\"id\": "<<(i+1)<<", \"center\": ["<<p.center[0]<<", "<<p.center[1]<<", "<<p.center[2]
                   <<"], \"volume_a3\": "<<p.volume<<", \"peak_buriedness\": "<<p.peakBuriedness
                   <<", \"points\": "<<p.pointCount
                   <<", \"mean_buriedness\": "<<p.meanBuriedness
                   <<", \"high_buriedness_fraction\": "<<p.highBuriednessFraction
                   <<", \"boundary_area_a2\": "<<p.exposedFaceArea
                   <<", \"compactness\": "<<p.compactness
                   <<", \"elongation\": "<<p.elongation
                   <<", \"mean_surface_distance_a\": "<<p.meanAtomSurfaceDistance
                   <<", \"heuristic_score\": "<<p.heuristicScore<<"}"<<(i+1<n?",":"")<<"\n";
        }
        *stream<<"  ]\n}\n";
        if(!points.empty()){
            pointDest.open(points);
            if(!pointDest)throw std::runtime_error("Cannot write: "+points);
            pointDest<<"pocket_id,x,y,z,buriedness\n"<<std::setprecision(8);
            for(std::size_t i=0;i<n;i++)for(const auto& p:r.pockets[i].points)
                pointDest<<i+1<<","<<p.pos[0]<<","<<p.pos[1]<<","<<p.pos[2]<<","<<p.buriedness<<"\n";
        }
        if(!pdb.empty()){
            std::size_t count=0;
            for(std::size_t i=0;i<n;i++){
                if(i>=9999)throw std::runtime_error("PDB supports at most 9999 pockets; use --top");
                if(r.pockets[i].points.size()>99999-count)
                    throw std::runtime_error("PDB exceeds 99999 pseudoatoms; increase --spacing or use --top");
                count+=r.pockets[i].points.size();
                for(const auto& pt:r.pockets[i].points)
                    for(double x:pt.pos)if(x<=-1000.0 || x>=10000.0)
                        throw std::runtime_error("Pocket coordinate is outside PDB 8.3 limits");
            }
            pdbDest.open(pdb);
            if(!pdbDest)throw std::runtime_error("Cannot write: "+pdb);
            pdbDest<<"REMARK 900 LIGINT POCKET BURIEDNESS PSEUDOATOMS\n"
                   <<"REMARK 900 B FACTOR = 100 * BURIEDNESS (0..100)\n"
                   <<"REMARK 900 RESNAME PKT CHAIN P RESSEQ = POCKET ID\n";
            std::size_t serial=0;
            for(std::size_t i=0;i<n;i++)for(const auto& pt:r.pockets[i].points){
                char line[160];
                const int written=std::snprintf(line,sizeof(line),
                    "HETATM%5u  C   PKT P%4u    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s",
                    static_cast<unsigned>(++serial),static_cast<unsigned>(i+1),
                    pt.pos[0],pt.pos[1],pt.pos[2],1.0,
                    100.0*static_cast<double>(pt.buriedness),"C");
                if(written<0 || written>=static_cast<int>(sizeof(line)))
                    throw std::runtime_error("Cannot format pocket PDB record");
                pdbDest<<line<<"\n";
            }
            pdbDest<<"END\n";
            if(!pdbDest)throw std::runtime_error("Failed writing PDB: "+pdb);
        }
    }catch(const std::exception& e){
        std::cerr<<"ligint-pocket: "<<e.what()<<"\n";return 1;
    }
    return 0;
}
