// Pocket input parser: uses the same PDB/mmCIF atom readers as the interaction analysis.
#include "internal.hpp"
#include "ligint/pocket.hpp"
#include <fstream>
#include <stdexcept>

namespace ligint {
std::vector<PocketAtom> readPocketAtoms(const std::string& path, StructureFormat fmt) {
    std::ifstream f(path, std::ios::binary);
    if(!f) throw std::runtime_error("Cannot open structure: "+path);
    if(fmt==StructureFormat::Auto) {
        // Infer from the first nonempty, non-comment line, including CIF data_ blocks.
        std::string line;
        while(std::getline(f,line)) {
            const auto a=line.find_first_not_of(" \t\r");
            if(a==std::string::npos || line[a]=='#')continue;
            fmt=(line.compare(a,5,"data_")==0 || line.compare(a,5,"loop_")==0 ||
                 line.compare(a,1,"_")==0) ? StructureFormat::Mmcif : StructureFormat::Pdb;
            break;
        }
        f.clear(); f.seekg(0);
    }
    std::vector<PocketAtom> out;
    auto sink=[&](detail::RawAtom& a) {
        // ATOM records only: excludes ligands, crystallographic water and most additives.
        if(a.het || a.el=="H" || a.el=="D")return;
        const auto el=detail::normEl(a.el);
        if(el=="H" || el=="D")return;
        for(double v:{a.x,a.y,a.z}) if(!std::isfinite(v))return;
        out.push_back({{a.x,a.y,a.z},detail::vdw(el),a.chain,
            a.res+" "+std::to_string(a.resi)+a.chain});
    };
    if(fmt==StructureFormat::Mmcif) {
        detail::CifStats stats;
        detail::readMmcifAtoms(f,sink,stats);
        if(!stats.sawAtomSite)throw std::runtime_error("No _atom_site table found in mmCIF");
    } else detail::readPdbAtoms(f,sink);
    if(out.empty())throw std::runtime_error("No protein ATOM records found");
    return out;
}
}
