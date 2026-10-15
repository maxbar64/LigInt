// ligint - command-line interface.
#include "ligint/fingerprint.hpp"
#include "ligint/ligint.hpp"
#include "ligint/klifs.hpp"
#include "ligint/rcsb.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace ligint;

static void usage(const char* prog, const Params& p) {
    cerr << "Usage: " << prog << " structure.pdb|structure.cif [RESNAME [chain [resseq[A|-]]]] [options]\n"
         << "Options: --sdf F  --mol2 F  --csv  --no-water  --no-merge  --plip-like  --metal-donors  --backbone-hyd  --altloc X  --format auto|pdb|cif  --version  --help\n"
         << "RCSB: --download-pdb ID [--download-format pdb|cif] [--download-output FILE] [--download-cache DIR] [--download-refresh] [--download-only]\n"
         << "Fingerprint: --fingerprint csv|json|svg  --pocket FILE  --klifs-pdb ID --klifs-chain C [--klifs-model M] [--klifs-id N] [--klifs-mapping FILE] [--klifs-output FILE]\n"
         << "Thresholds (defaults): --hb-min " << p.hbMin << " --hb-max " << p.hbMax << " --hb-ha " << p.hbHA
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
    for (int i = 1; i < argc; i++)            // preset applied first: explicit options override it
        if (string(argv[i]) == "--plip-like") P.applyPlipLike();
    vector<string> pa;
    bool csv = false;
    string fpFormat, pocketPath, klifsOutput; KlifsOptions klifs; bool useKlifs=false; RcsbDownloadOptions rcsb; bool useRcsb=false, downloadOnly=false;
    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "--help" || s == "-h") { usage(argv[0], P); return 0; }
        else if (s == "--version") { cout << "ligint " << LIGINT_VERSION << "\n"; return 0; }
        else if (s == "--csv") csv = true;
        else if (s == "--no-water") opt.useWater = false;
        else if (s == "--metal-donors") opt.metalDonors = true;
        else if (s == "--backbone-hyd") opt.backboneHyd = true;
        else if (s == "--no-merge") opt.mergeLinked = false;
        else if (s == "--plip-like") {}           // already applied above
        else if (s == "--format") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            const string v = argv[++i];
            if (v == "auto") opt.format = StructureFormat::Auto;
            else if (v == "pdb") opt.format = StructureFormat::Pdb;
            else if (v == "cif" || v == "mmcif") opt.format = StructureFormat::Mmcif;
            else { cerr << "--format expects auto, pdb or cif\n"; return 1; }
        }
        else if (s == "--altloc") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            const string v = argv[++i];
            if (v.size() != 1 || v[0] == ' ') { cerr << "--altloc expects a single character, e.g. A or B\n"; return 1; }
            opt.altloc = v[0];
        }
        else if (s == "--klifs-pdb" || s == "--klifs-chain" || s == "--klifs-model" || s == "--klifs-mapping" || s == "--klifs-output" || s == "--klifs-id") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            string v=argv[++i];useKlifs=true;
            if(s=="--klifs-pdb")klifs.pdbId=v;
            else if(s=="--klifs-chain")klifs.chain=v;
            else if(s=="--klifs-model")klifs.model=v;
            else if(s=="--klifs-mapping")klifs.mappingFile=v;
            else if(s=="--klifs-output")klifsOutput=v;
            else {try{size_t n=0;klifs.structureId=stoi(v,&n);if(n!=v.size()||klifs.structureId<=0)throw invalid_argument("id");}catch(...){cerr<<"Invalid KLIFS structure ID\n";return 1;}}
        }
        else if(s=="--download-refresh") {rcsb.refresh=true;}
        else if(s=="--download-only") {downloadOnly=true;}
        else if(s=="--download-pdb" || s=="--download-format" || s=="--download-output" || s=="--download-cache") {
            if(i+1>=argc){cerr<<s<<" requires a value\n";return 1;}
            string v=argv[++i];
            if(s=="--download-pdb"){rcsb.pdbId=v;useRcsb=true;}
            else if(s=="--download-format") rcsb.format=v;
            else if(s=="--download-output") rcsb.outputFile=v;
            else rcsb.cacheDir=v;
        }
        else if (s == "--fingerprint" || s == "--pocket") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            string v = argv[++i];
            if (s == "--pocket") pocketPath = v;
            else {
                if (v != "csv" && v != "json" && v != "svg") { cerr << "--fingerprint expects csv, json or svg\n"; return 1; }
                fpFormat = v;
            }
        }
        else if (s == "--sdf" || s == "--mol2") {
            if (i + 1 >= argc) { cerr << s << " requires a file\n"; return 1; }
            if (!ifstream(argv[i + 1])) { cerr << "Cannot open " << argv[i + 1] << "\n"; return 1; }
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
                cerr << "Unknown option, missing value or non-numeric: " << s << "\n"; usage(argv[0], P); return 1;
            }
            i++;
        } else pa.push_back(s);
    }
    if(downloadOnly && !useRcsb){cerr<<"--download-only requires --download-pdb\n";return 1;}
    // with --download-pdb the positional arguments are the ligand: RESNAME [chain [resseq]] (the structure is the downloaded file)
    if(useRcsb && !pa.empty()){ifstream existing(pa[0]);if(existing){cerr<<"Use --download-pdb without a separate structure path (the arguments after the options are the ligand: RESNAME [chain [resseq]])\n";return 1;}}
    if(useRcsb){try{string path=downloadRcsbStructure(rcsb);cerr<<"RCSB structure: "<<path<<"\n";if(downloadOnly){cout<<path<<"\n";return 0;}pa.insert(pa.begin(),path);}catch(const exception& e){cerr<<"RCSB: "<<e.what()<<"\n";return 1;}}
    if (pa.empty() && !useKlifs) { usage(argv[0], P); return 1; }
    for (auto& w : P.clampToGrid()) cerr << w << "\n";
    opt.ligName = pa.size() > 1 ? pa[1] : "";
    opt.chain = pa.size() > 2 ? pa[2] : "";
    if (pa.size() > 3) {
        // residue number with an optional insertion code: "501" (any code), "501A" (code A), "501-" (no code)
        const string& rs = pa[3];
        size_t k = 0;
        if (k < rs.size() && (rs[k] == '-' || rs[k] == '+')) k++;
        const size_t d0 = k;
        while (k < rs.size() && isdigit(static_cast<unsigned char>(rs[k]))) k++;
        const string suffix = rs.substr(k);
        const bool okSuffix = suffix.empty() || suffix == "-" || (suffix.size() == 1 && isalpha(static_cast<unsigned char>(suffix[0])));
        if (k == d0 || k - d0 > 6 || !okSuffix) { cerr << "Invalid residue number: " << rs << " (expected e.g. 501, 501A or 501-)\n"; return 1; }
        opt.resi = atoi(rs.substr(0, k).c_str());
        opt.icode = suffix.empty() ? 0 : (suffix == "-" ? ' ' : suffix[0]);
    }

    PocketDefinition pocket;
    if(useKlifs && !pocketPath.empty()){cerr<<"Choose either --pocket or --klifs-*\n";return 1;}
    if(useKlifs && pa.empty() && klifsOutput.empty()){cerr<<"--klifs-* without a structure needs --klifs-output FILE (otherwise the pocket would not be saved anywhere)\n";return 1;}
    if(useKlifs && !pa.empty() && fpFormat.empty()) cerr<<"Warning: the KLIFS pocket is only used with --fingerprint\n";
    if(useKlifs){try{pocket=generateKlifsPocket(klifs);if(!klifsOutput.empty())saveKlifsPocket(pocket,klifsOutput);cerr<<"KLIFS pocket: "<<pocket.positions.size()<<" mapped positions\n";}catch(const exception& e){cerr<<"KLIFS: "<<e.what()<<"\n";return 1;}}
    if(pa.empty()) return 0; // generate TSV without structure analysis
    const bool havePocket = useKlifs || !pocketPath.empty();
    if (!pocketPath.empty()) {
        vector<string> pw;
        if (!loadPocketDefinition(pocketPath, pocket, pw)) { cerr << "Cannot open " << pocketPath << "\n"; return 1; }
        for (auto& w : pw) cerr << "Warning: " << w << "\n";
    }
    if (havePocket && fpFormat.empty() && !useKlifs) { cerr << "--pocket requires --fingerprint\n"; return 1; }

    try {
        Result r = analyzeFile(pa[0], opt);
        for (auto& w : r.warnings) cerr << w << "\n";
        if (!fpFormat.empty() && !r.listingOnly) {
            vector<Fingerprint> fps;
            for (auto& lig : r.ligands) {
                fps.push_back(makeFingerprint(lig, havePocket ? &pocket : nullptr));
                for (auto& w : fps.back().warnings) cerr << "Warning: " << w << "\n";
            }
            for (auto& s : r.info) cerr << s << "\n";
            for (auto& lig : r.ligands) for (auto& s : lig.info) cerr << s << "\n";
            if (fpFormat == "csv") writeFingerprintCsv(cout, fps);
            else if (fpFormat == "json") writeFingerprintJson(cout, fps);
            else {
                if (fps.size() > 1) cerr << "Warning: " << fps.size() << " ligands found, only the first is drawn (select one with chain and residue)\n";
                cout << fingerprintSvg(fps.front());
            }
        } else if (csv) printCsv(cout, cerr, r);
        else printText(cout, r);
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
