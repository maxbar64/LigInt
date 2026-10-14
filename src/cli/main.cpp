// ligint - interfaccia a riga di comando.
#include "ligint/ligint.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace ligint;

static void usage(const char* prog, const Params& p) {
    cerr << "Uso: " << prog << " file.pdb [RESNAME [catena [resseq]]] [opzioni]\n"
         << "Opzioni: --sdf F  --mol2 F  --csv  --no-water  --no-merge  --plip-like  --metal-donors  --backbone-hyd  --version  --help\n"
         << "Soglie (default): --hb-min " << p.hbMin << " --hb-max " << p.hbMax << " --hb-ha " << p.hbHA
         << " --hb-angle " << p.hbDHA << "\n  --hyd-min " << p.hydMin << " --hyd-max " << p.hydMax
         << " --salt-max " << p.saltMax << " --rep-max " << p.repMax << "\n  --hal-max " << p.halMax
         << " --hal-angle " << p.halAngle << " --metal-max " << p.metalMax << " --pi-max " << p.piMax
         << " --pi-offset " << p.piOffset << " --pi-angle " << p.piAngle << " --cp-max " << p.cpMax << " --clash " << p.clash
         << "\n  --wb-min " << p.wbMin << " --wb-max " << p.wbMax << " --hb-ante " << p.hbAnte
         << " --hb-long " << p.hbLong << " --hb-long-tol " << p.hbLongTol << " --salt-centroid " << p.saltCentroid << "\n";
}

int main(int argc, char** argv) {
    Options opt;
    Params& P = opt.params;
    for (int i = 1; i < argc; i++)            // preset applicato per primo: le opzioni esplicite lo sovrascrivono
        if (string(argv[i]) == "--plip-like") P.applyPlipLike();
    vector<string> pa;
    bool csv = false;
    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "--help" || s == "-h") { usage(argv[0], P); return 0; }
        else if (s == "--version") { cout << "ligint " << LIGINT_VERSION << "\n"; return 0; }
        else if (s == "--csv") csv = true;
        else if (s == "--no-water") opt.useWater = false;
        else if (s == "--metal-donors") opt.metalDonors = true;
        else if (s == "--backbone-hyd") opt.backboneHyd = true;
        else if (s == "--no-merge") opt.mergeLinked = false;
        else if (s == "--plip-like") {}           // gia' applicato sopra
        else if (s == "--sdf" || s == "--mol2") {
            if (i + 1 >= argc) { cerr << s << " richiede un file\n"; return 1; }
            if (!ifstream(argv[i + 1])) { cerr << "Impossibile aprire " << argv[i + 1] << "\n"; return 1; }
            vector<string> warnings;
            bool ok = s == "--sdf" ? loadSdf(argv[i + 1], opt.templates, warnings) : loadMol2(argv[i + 1], opt.templates, warnings);
            i++;
            for (auto& w : warnings) cerr << w << "\n";
            if (!ok) return 1;
        } else if (s.compare(0, 2, "--") == 0) {
            char* endp = nullptr;
            double val = i + 1 < argc ? strtod(argv[i + 1], &endp) : 0.0;
            bool numOK = i + 1 < argc && endp != argv[i + 1] && *endp == '\0' && isfinite(val);
            if (!numOK || !P.set(s.substr(2), val)) {
                cerr << "Opzione sconosciuta, senza valore o non numerica: " << s << "\n"; usage(argv[0], P); return 1;
            }
            i++;
        } else pa.push_back(s);
    }
    if (pa.empty()) { usage(argv[0], P); return 1; }
    for (auto& w : P.clampToGrid()) cerr << w << "\n";
    opt.ligName = pa.size() > 1 ? pa[1] : "";
    opt.chain = pa.size() > 2 ? pa[2] : "";
    opt.resi = pa.size() > 3 ? atoi(pa[3].c_str()) : INT_MIN;

    try {
        Result r = analyzeFile(pa[0], opt);
        for (auto& w : r.warnings) cerr << w << "\n";
        if (csv) printCsv(cout, cerr, r);
        else printText(cout, r);
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
