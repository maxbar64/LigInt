// Reader of the legacy PDB format: the ATOM and HETATM records of the first model.
#include "internal.hpp"

namespace ligint {
namespace detail {

void readPdbAtoms(std::istream& in, const AtomSink& sink) {
    using namespace std;
    string line;
    while (getline(in, line)) {
        if (line.compare(0, 6, "ENDMDL") == 0) break;  // first model only
        if (line.size() < 54) continue;
        const string rec = line.substr(0, 6);
        if (rec != "ATOM  " && rec != "HETATM") continue;
        RawAtom r;
        r.name = trim(line.substr(12, 4));
        r.res = trim(line.substr(17, 3));
        r.chain = string(1, line[21]);
        r.resi = atoi(line.substr(22, 4).c_str());
        r.icode = line[26];
        r.alt = line[16];
        r.x = atof(line.substr(30, 8).c_str());
        r.y = atof(line.substr(38, 8).c_str());
        r.z = atof(line.substr(46, 8).c_str());
        r.bf = line.size() >= 66 ? atof(line.substr(60, 6).c_str()) : 0.0;
        r.het = rec == "HETATM";
        string el = line.size() >= 78 ? trim(line.substr(76, 2)) : "";
        if (el.empty()) {
            const string nm = line.substr(12, 4);
            if (r.het && r.name.size() <= 2 && isalpha((unsigned char)nm[0]) && isalpha((unsigned char)nm[1]))
                el = nm.substr(0, 2);
            else for (char c : r.name) if (isalpha((unsigned char)c)) { el = string(1, c); break; }
        }
        r.el = el;
        sink(r);
    }
}

}  // namespace detail
}  // namespace ligint
