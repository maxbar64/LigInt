// ligint - command-line interface.
#include "ligint/fingerprint.hpp"
#include "ligint/ligint.hpp"
#include "ligint/gpcrdb.hpp"
#include "ligint/klifs.hpp"
#include "ligint/matrix.hpp"
#include "ligint/rcsb.hpp"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

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
         << "Options: --sdf F  --mol2 F  --exclude-pdb-ligand RESNAME  --external-index N (1-based)  --csv  --no-water  --no-merge  --plip-like  --metal-donors  --backbone-hyd  --altloc X  --format auto|pdb|cif  --version  --help\n"
         << "RCSB: --download-pdb ID [--download-format pdb|cif] [--download-output FILE] [--download-cache DIR] [--download-refresh] [--download-only]\n"
         << "Fingerprint: --fingerprint csv|json|svg|OUTPUT.csv|OUTPUT.json|OUTPUT.svg  --pocket FILE  --klifs-pdb ID --klifs-chain C [--klifs-model M] [--klifs-id N] [--klifs-mapping FILE] [--klifs-sequence S] [--klifs-output FILE]\n"
         << "GPCR pocket: --gpcrdb-positions LIST|FILE [--gpcrdb-pdb FILE | --gpcrdb-mapping FILE] [--gpcrdb-chain C] [--gpcrdb-name TEXT] [--gpcrdb-output FILE]\n"
         << "Matrix: --matrix LIST [--matrix-format csv|json|svg] [--matrix-classes klifs|all|NAMES] [--matrix-hide-empty] [--matrix-title TEXT] [--pocket FILE]\n"
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
    string matrixPath, matrixFormat = "csv", matrixClasses, matrixTitle; bool matrixHide = false, matrixOptionGiven = false;
    string fpFormat, fpOutput, pocketPath, klifsOutput; KlifsOptions klifs; bool useKlifs=false; GpcrdbOptions gpcrdb; bool useGpcrdb=false; string gpcrdbOutput, gpcrdbPositions; RcsbDownloadOptions rcsb; bool useRcsb=false, downloadOnly=false;
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
        else if (s == "--matrix" || s == "--matrix-format" || s == "--matrix-classes" || s == "--matrix-title") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            const string v = argv[++i];
            if (s == "--matrix") matrixPath = v;
            else {
                matrixOptionGiven = true;
                if (s == "--matrix-format") {
                    if (v != "csv" && v != "json" && v != "svg") { cerr << "--matrix-format expects csv, json or svg\n"; return 1; }
                    matrixFormat = v;
                } else if (s == "--matrix-classes") matrixClasses = v;
                else matrixTitle = v;
            }
        }
        else if (s == "--matrix-hide-empty") { matrixHide = true; matrixOptionGiven = true; }
        else if (s == "--gpcrdb-positions" || s == "--gpcrdb-pdb" || s == "--gpcrdb-mapping" || s == "--gpcrdb-chain" || s == "--gpcrdb-name" || s == "--gpcrdb-output") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            const string v = argv[++i];
            useGpcrdb = true;
            if (s == "--gpcrdb-positions") gpcrdbPositions = v;
            else if (s == "--gpcrdb-pdb") gpcrdb.numberedStructure = v;
            else if (s == "--gpcrdb-mapping") gpcrdb.mappingFile = v;
            else if (s == "--gpcrdb-chain") gpcrdb.chain = v;
            else if (s == "--gpcrdb-name") gpcrdb.name = v;
            else gpcrdbOutput = v;
        }
        else if (s == "--klifs-pdb" || s == "--klifs-chain" || s == "--klifs-model" || s == "--klifs-mapping" || s == "--klifs-output" || s == "--klifs-id" || s == "--klifs-sequence") {
            if (i + 1 >= argc) { cerr << s << " requires a value\n"; return 1; }
            string v=argv[++i];useKlifs=true;
            if(s=="--klifs-pdb")klifs.pdbId=v;
            else if(s=="--klifs-chain")klifs.chain=v;
            else if(s=="--klifs-model")klifs.model=v;
            else if(s=="--klifs-mapping")klifs.mappingFile=v;
            else if(s=="--klifs-output")klifsOutput=v;
            else if(s=="--klifs-sequence")klifs.pocketSequence=v;
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
                if (v == "csv" || v == "json" || v == "svg") fpFormat = v;
                else {
                    const auto dot = v.find_last_of('.');
                    const string ext = dot == string::npos ? "" : v.substr(dot + 1);
                    if (ext != "csv" && ext != "json" && ext != "svg") {
                        cerr << "--fingerprint expects csv, json, svg, or a file ending in .csv, .json or .svg\n";
                        return 1;
                    }
                    fpFormat = ext;
                    fpOutput = v;
                }
            }
        }
        else if (s == "--external-index") {
            if (i + 1 >= argc) { cerr << "--external-index requires a positive integer\n"; return 1; }
            const std::string value = argv[++i];
            char* end = nullptr;
            const long index = std::strtol(value.c_str(), &end, 10);
            if (value.empty() || *end != '\0' || index < 1 || index > 1000000) {
                cerr << "--external-index requires a positive integer (1-based)\n"; return 1;
            }
            opt.externalLigandIndex = static_cast<int>(index - 1);
        }
        else if (s == "--exclude-pdb-ligand") {
            if (i + 1 >= argc || argv[i + 1][0] == '-') { cerr << "--exclude-pdb-ligand requires a residue name\n"; return 1; }
            opt.excludePdbLigand = argv[++i];
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
    if (matrixPath.empty() && matrixOptionGiven) { cerr << "--matrix-format, --matrix-classes, --matrix-title and --matrix-hide-empty require --matrix LIST\n"; return 1; }
    if (!matrixPath.empty()) {
        // the matrix: one row for each ligand of the list, the structures are in the list
        if (!pa.empty()) { cerr << "--matrix takes no structure argument: the structures are in the list\n"; return 1; }
        if (useRcsb || useKlifs || useGpcrdb || !fpFormat.empty() || csv) { cerr << "--matrix cannot be combined with --download-pdb, --klifs-*, --gpcrdb-*, --fingerprint or --csv\n"; return 1; }
        for (auto& w : P.clampToGrid()) cerr << w << "\n";
        try {
            MatrixOptions mo;
            mo.analysis = opt;
            mo.pocket = pocketPath;                       // the pocket of the entries that do not name one
            mo.hideEmptyPositions = matrixHide;
            if (!matrixClasses.empty()) {
                string err;
                if (!parseMatrixClasses(matrixClasses, mo.classes, err)) { cerr << "--matrix-classes: " << err << "\n"; return 1; }
            }
            const InteractionMatrix m = makeMatrix(readMatrixList(matrixPath), mo);
            for (auto& w : m.warnings) cerr << "Warning: " << w << "\n";
            if (matrixFormat == "csv") writeMatrixCsv(cout, m);
            else if (matrixFormat == "json") writeMatrixJson(cout, m);
            else { MatrixSvgOptions so; so.title = matrixTitle; cout << matrixSvg(m, so); }
        } catch (const exception& e) {
            cerr << "Matrix: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if(downloadOnly && !useRcsb){cerr<<"--download-only requires --download-pdb\n";return 1;}
    // with --download-pdb the positional arguments are the ligand: RESNAME [chain [resseq]] (the structure is the downloaded file)
    if(useRcsb && !pa.empty()){ifstream existing(pa[0]);if(existing){cerr<<"Use --download-pdb without a separate structure path (the arguments after the options are the ligand: RESNAME [chain [resseq]])\n";return 1;}}
    if(useRcsb){try{string path=downloadRcsbStructure(rcsb);cerr<<"RCSB structure: "<<path<<"\n";if(downloadOnly){cout<<path<<"\n";return 0;}pa.insert(pa.begin(),path);}catch(const exception& e){cerr<<"RCSB: "<<e.what()<<"\n";return 1;}}
    if (pa.empty() && !useKlifs && !useGpcrdb) { usage(argv[0], P); return 1; }
    for (auto& w : P.clampToGrid()) cerr << w << "\n";
    opt.ligName = pa.size() > 1 ? pa[1] : "";
    if (pa.size() == 1 && !opt.templates.empty() && opt.externalLigandIndex < 0) opt.externalLigandIndex = 0;
    if (opt.externalLigandIndex >= 0 && (pa.size() != 1 || static_cast<size_t>(opt.externalLigandIndex) >= opt.templates.size())) {
        cerr << "--external-index needs a receptor without a RESNAME argument and a valid loaded SDF/MOL2 molecule index\n"; return 1;
    }
    if (!opt.excludePdbLigand.empty() && opt.externalLigandIndex < 0) {
        cerr << "--exclude-pdb-ligand requires an external SDF/MOL2 ligand pose (omit RESNAME)\n"; return 1;
    }
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
    if(useGpcrdb && useKlifs){cerr<<"Choose either --klifs-* or --gpcrdb-*\n";return 1;}
    if(useGpcrdb && !pocketPath.empty()){cerr<<"Choose either --pocket or --gpcrdb-*\n";return 1;}
    if(useKlifs && !pocketPath.empty()){cerr<<"Choose either --pocket or --klifs-*\n";return 1;}
    if(useKlifs && pa.empty() && klifsOutput.empty()){cerr<<"--klifs-* without a structure needs --klifs-output FILE (otherwise the pocket would not be saved anywhere)\n";return 1;}
    if(useKlifs && !pa.empty() && fpFormat.empty()) cerr<<"Warning: the KLIFS pocket is only used with --fingerprint\n";
    if(useKlifs){try{pocket=generateKlifsPocket(klifs);if(!klifsOutput.empty())saveKlifsPocket(pocket,klifsOutput);cerr<<"KLIFS pocket: "<<pocket.positions.size()<<" mapped positions\n";}catch(const exception& e){cerr<<"KLIFS: "<<e.what()<<"\n";return 1;}}
    if(useGpcrdb && pa.empty() && gpcrdbOutput.empty()){cerr<<"--gpcrdb-* without a structure needs --gpcrdb-output FILE (otherwise the pocket would not be saved anywhere)\n";return 1;}
    if(useGpcrdb && !pa.empty() && fpFormat.empty()) cerr<<"Warning: the GPCRdb pocket is only used with --fingerprint\n";
    if(useGpcrdb){
        try{
            if(gpcrdbPositions.empty()){cerr<<"--gpcrdb-positions is required (the pocket: generic residue numbers, e.g. 3x32,6x48 or a file with them)\n";return 1;}
            gpcrdb.positions=parseGpcrdbPositions(gpcrdbPositions);
            if(gpcrdb.numberedStructure.empty() && gpcrdb.mappingFile.empty()){
                // the structure of the analysis can be the file annotated by GPCRdb
                if(pa.empty()){cerr<<"--gpcrdb-* without a structure needs --gpcrdb-pdb FILE (annotated by GPCRdb) or --gpcrdb-mapping FILE\n";return 1;}
                gpcrdb.numberedStructure=pa[0];
            }
            const GpcrdbPocket gp=generateGpcrdbPocket(gpcrdb);
            for(auto& w:gp.warnings)cerr<<"Warning: "<<w<<"\n";
            pocket=gp.pocket;
            if(!gpcrdbOutput.empty())saveGpcrdbPocket(pocket,gpcrdbOutput);
            cerr<<"GPCRdb pocket: "<<pocket.positions.size()<<" positions\n";
        }catch(const exception& e){cerr<<"GPCRdb: "<<e.what()<<"\n";return 1;}
    }
    if(pa.empty()) return 0; // generate TSV without structure analysis
    const bool havePocket = useKlifs || useGpcrdb || !pocketPath.empty();
    if (!pocketPath.empty()) {
        vector<string> pw;
        if (!loadPocketDefinition(pocketPath, pocket, pw)) { cerr << "Cannot open " << pocketPath << "\n"; return 1; }
        for (auto& w : pw) cerr << "Warning: " << w << "\n";
    }
    if (havePocket && fpFormat.empty() && !useKlifs && !useGpcrdb) { cerr << "--pocket requires --fingerprint\n"; return 1; }

    try {
        Result r = analyzeFile(pa[0], opt);
        for (auto& w : r.warnings) cerr << w << "\n";
        #ifdef _WIN32
        if (fpOutput.empty()) {
            if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
                cerr << "Cannot set binary stdout mode\n";
                return 1;
            }
        }
        #endif
        if (!fpFormat.empty() && !r.listingOnly) {
            vector<Fingerprint> fps;
            for (auto& lig : r.ligands) {
                fps.push_back(makeFingerprint(lig, havePocket ? &pocket : nullptr));
                for (auto& w : fps.back().warnings) cerr << "Warning: " << w << "\n";
            }
            for (auto& s : r.info) cerr << s << "\n";
            for (auto& lig : r.ligands) for (auto& s : lig.info) cerr << s << "\n";
            ofstream outputFile;
            ostream* destination = &cout;
            if (!fpOutput.empty()) {
                outputFile.open(fpOutput, ios::binary | ios::trunc);
                if (!outputFile) { cerr << "Cannot open fingerprint output file: " << fpOutput << "\n"; return 1; }
                destination = &outputFile;
            }
            if (fpFormat == "csv") writeFingerprintCsv(*destination, fps);
            else if (fpFormat == "json") writeFingerprintJson(*destination, fps);
            else {
                if (fps.size() > 1) cerr << "Warning: " << fps.size() << " ligands found, only the first is drawn (select one with chain and residue)\n";
                if (fps.empty()) { cerr << "No ligand fingerprint available for SVG output\n"; return 1; }
                *destination << fingerprintSvg(fps.front());
            }
            if (!*destination) { cerr << "Failed to write fingerprint output" << (fpOutput.empty() ? "" : ": " + fpOutput) << "\n"; return 1; }
        } else if (csv) printCsv(cout, cerr, r);
        else printText(cout, r);
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }
    return 0;
}
