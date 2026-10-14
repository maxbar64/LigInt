// Interaction fingerprints: pocket definition, per-residue classification, CSV and JSON output.
#include "ligint/fingerprint.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>

namespace ligint {

// ---------------------------------------------------------------------------------------------
// classes
// ---------------------------------------------------------------------------------------------
const char* fpClassName(int c) {
    static const char* n[FpClassCount] = {"hydrophobic", "face_to_face", "edge_to_face", "hbond_donor", "hbond_acceptor",
                                          "ionic_positive", "ionic_negative", "cation_pi", "halogen", "metal", "water_bridge"};
    return c >= 0 && c < FpClassCount ? n[c] : "?";
}

const char* fpClassTitle(int c) {
    static const char* n[FpClassCount] = {"Hydrophobic", "Aromatic face-to-face", "Aromatic edge-to-face", "H-bond donor",
                                          "H-bond acceptor", "Ionic positive", "Ionic negative", "Cation-pi", "Halogen bond",
                                          "Metal coordination", "Water bridge"};
    return c >= 0 && c < FpClassCount ? n[c] : "?";
}

bool FingerprintEntry::any() const {
    for (int c : counts) if (c > 0) return true;
    return false;
}

char oneLetterCode(const std::string& r) {
    static const std::map<std::string, char> m = {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
        {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'}, {"SER", 'S'},
        {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}, {"MSE", 'M'}, {"SEP", 'S'}, {"TPO", 'T'}, {"PTR", 'Y'},
        {"DA", 'A'}, {"DC", 'C'}, {"DG", 'G'}, {"DT", 'T'}, {"A", 'A'}, {"C", 'C'}, {"G", 'G'}, {"U", 'U'}};
    auto it = m.find(r);
    return it == m.end() ? '?' : it->second;
}

// ---------------------------------------------------------------------------------------------
// pocket definition
// ---------------------------------------------------------------------------------------------
namespace {

std::string trimStr(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}

// splits on tabs or commas when present (so that fields may contain spaces), otherwise on whitespace
std::vector<std::string> splitFields(const std::string& line) {
    std::vector<std::string> f;
    char sep = line.find('\t') != std::string::npos ? '\t' : (line.find(',') != std::string::npos ? ',' : ' ');
    if (sep == ' ') {
        std::istringstream ss(line);
        std::string w;
        while (ss >> w) f.push_back(w);
    } else {
        std::string cur;
        for (char c : line) {
            if (c == sep) { f.push_back(trimStr(cur)); cur.clear(); }
            else cur += c;
        }
        f.push_back(trimStr(cur));
    }
    return f;
}

// "134", "134A", "-5" -> number and insertion code; false if malformed
bool parseResi(const std::string& s, int& resi, char& icode) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[i] == '-' || s[i] == '+') i++;
    size_t d0 = i;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) i++;
    if (i == d0 || i - d0 > 6) return false;
    resi = std::atoi(s.substr(0, i).c_str());
    icode = ' ';
    if (i < s.size()) {
        if (i + 1 != s.size() || !std::isalpha(static_cast<unsigned char>(s[i]))) return false;
        icode = s[i];
    }
    return true;
}

}  // namespace

bool parsePocketDefinition(std::istream& in, PocketDefinition& out, std::vector<std::string>& warnings) {
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        n++;
        std::string t = trimStr(line);
        if (t.empty() || t[0] == '#') continue;
        if (t[0] == '@') {
            std::istringstream ss(t);
            std::string key;
            ss >> key;
            std::string rest;
            std::getline(ss, rest);
            rest = trimStr(rest);
            if (key == "@name") out.name = rest;
            else if (key == "@color") {
                std::istringstream cs(rest);
                std::string region, color;
                cs >> region >> color;
                if (region.empty() || color.size() != 7 || color[0] != '#') warnings.push_back("pocket definition line " + std::to_string(n) + ": bad @color");
                else out.colors.push_back({region, color});
            } else warnings.push_back("pocket definition line " + std::to_string(n) + ": unknown directive " + key);
            continue;
        }
        std::vector<std::string> f = splitFields(t);
        // chain, residue number, position label, region [, residue name]
        if (f.size() < 4) { warnings.push_back("pocket definition line " + std::to_string(n) + ": expected chain, residue, position, region"); continue; }
        PocketPosition p;
        if (!parseResi(f[1], p.resi, p.icode)) { warnings.push_back("pocket definition line " + std::to_string(n) + ": bad residue number '" + f[1] + "'"); continue; }
        p.chain = (f[0] == "*" || f[0] == "-" || f[0] == "_") ? "" : f[0];
        p.label = f[2];
        p.region = f[3];
        if (f.size() >= 5) {
            p.resName = f[4];
            std::transform(p.resName.begin(), p.resName.end(), p.resName.begin(), [](unsigned char c) { return std::toupper(c); });
        }
        out.positions.push_back(p);
    }
    return true;
}

bool loadPocketDefinition(const std::string& path, PocketDefinition& out, std::vector<std::string>& warnings) {
    std::ifstream f(path);
    if (!f) return false;
    return parsePocketDefinition(f, out, warnings);
}

// ---------------------------------------------------------------------------------------------
// fingerprint
// ---------------------------------------------------------------------------------------------
namespace {

struct ResId {
    std::string resName, chain;
    int resi = 0;
    char icode = ' ';
};

// "ASP 189A" / "ASP 189AA" (insertion code + chain) -> parts. Labels are built by the analysis as
// resName + " " + number + [icode] + chain.
bool parseResidueLabel(const std::string& label, ResId& id) {
    size_t sp = label.find(' ');
    if (sp == std::string::npos || sp == 0) return false;
    id.resName = label.substr(0, sp);
    size_t i = sp + 1, d0 = i;
    if (i < label.size() && label[i] == '-') i++;
    while (i < label.size() && std::isdigit(static_cast<unsigned char>(label[i]))) i++;
    if (i == d0) return false;
    id.resi = std::atoi(label.substr(d0, i - d0).c_str());
    std::string rest = label.substr(i);          // [icode] + chain (chain may be a blank)
    if (rest.empty()) { id.chain = ""; id.icode = ' '; }
    else if (rest.size() == 1) { id.chain = rest; id.icode = ' '; }
    else { id.icode = rest[0]; id.chain = rest.substr(1); }
    if (id.chain == " ") id.chain = "";
    return true;
}

bool sameResidue(const FingerprintEntry& e, const ResId& id) {
    return e.resi == id.resi && e.icode == id.icode && e.chain == id.chain;
}

void addCount(FingerprintEntry& e, int cls, double d) {
    e.counts[cls]++;
    if (e.minDistance[cls] == 0 || d < e.minDistance[cls]) e.minDistance[cls] = d;
}

// classes of one interaction (an interaction can feed more than one class, e.g. a charge-assisted H-bond)
std::vector<int> classesOf(const HitRecord& h) {
    std::vector<int> c;
    if (!h.favorable) return c;
    const std::string& t = h.type;
    if (t == type::Hydrophobic || t == type::BackboneApolar) c.push_back(FpHydrophobic);
    else if (t == type::PiStackingParallel) c.push_back(FpFaceToFace);
    else if (t == type::PiStackingT) c.push_back(FpEdgeToFace);
    else if (t == type::CationPi) c.push_back(FpCationPi);
    else if (t == type::HalogenBond) c.push_back(FpHalogen);
    else if (t == type::MetalCoordination) c.push_back(FpMetal);
    else if (t == type::WaterBridge) c.push_back(FpWaterBridge);
    else if (t == type::HydrogenBond || t == type::SaltBridge) {
        if (h.hbResidueRole == 'D') c.push_back(FpHBondDonor);
        else if (h.hbResidueRole == 'A') c.push_back(FpHBondAcceptor);
        if (t == type::SaltBridge) {
            if (h.residueCharge > 0) c.push_back(FpIonicPositive);
            else if (h.residueCharge < 0) c.push_back(FpIonicNegative);
        }
    }
    return c;
}

}  // namespace

Fingerprint makeFingerprint(const LigandResult& lig, const PocketDefinition* pocket) {
    Fingerprint fp;
    fp.ligandId = lig.id;
    if (pocket) { fp.pocketName = pocket->name; fp.regionColors = pocket->colors; }

    // pocket positions first (in definition order); contact residues are then matched or appended
    std::vector<int> posIndex;                 // entries[i] <- position index (or -1 for extra residues)
    if (pocket) {
        for (const auto& p : pocket->positions) {
            FingerprintEntry e;
            e.chain = p.chain; e.resi = p.resi; e.icode = p.icode;
            e.resName = p.resName;
            e.position = p.label; e.region = p.region;
            e.residue = (p.resName.empty() ? "?" : p.resName) + " " + std::to_string(p.resi) + (p.icode != ' ' ? std::string(1, p.icode) : "") + p.chain;
            fp.entries.push_back(e);
        }
    }
    const size_t nPocket = fp.entries.size();

    for (const auto& h : lig.hits) {
        std::vector<int> cls = classesOf(h);
        if (cls.empty()) continue;
        ResId id;
        if (!parseResidueLabel(h.residue, id)) continue;
        FingerprintEntry* target = nullptr;
        for (size_t i = 0; i < nPocket && !target; i++) {
            FingerprintEntry& e = fp.entries[i];
            const bool chainOk = e.chain.empty() || e.chain == id.chain;
            if (chainOk && e.resi == id.resi && e.icode == id.icode) target = &e;
        }
        if (target) {
            if (target->chain.empty()) target->chain = id.chain;
            target->residue = h.residue;
            if (target->resName.empty()) {
                target->resName = id.resName;
            } else if (target->resName != id.resName) {
                const std::string w = "residue " + h.residue + " does not match the pocket definition (" + target->resName + " expected)";
                if (std::find(fp.warnings.begin(), fp.warnings.end(), w) == fp.warnings.end()) fp.warnings.push_back(w);
            }
        } else {
            for (size_t i = nPocket; i < fp.entries.size() && !target; i++)
                if (sameResidue(fp.entries[i], id)) target = &fp.entries[i];
            if (!target) {
                FingerprintEntry e;
                e.residue = h.residue; e.resName = id.resName; e.chain = id.chain; e.resi = id.resi; e.icode = id.icode;
                e.inPocket = pocket == nullptr;
                fp.entries.push_back(e);
                target = &fp.entries.back();
            }
        }
        for (int c : cls) addCount(*target, c, h.distance);
    }

    // without a definition: sorted by chain, number, insertion code; with one: extra residues sorted after the pocket
    auto key = [](const FingerprintEntry& e) { return std::make_tuple(e.chain, e.resi, e.icode); };
    std::stable_sort(fp.entries.begin() + nPocket, fp.entries.end(),
                     [&](const FingerprintEntry& a, const FingerprintEntry& b) { return key(a) < key(b); });
    return fp;
}

// ---------------------------------------------------------------------------------------------
// CSV and JSON
// ---------------------------------------------------------------------------------------------
namespace {

std::string csvField(const std::string& s) {
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string r = "\"";
    for (char c : s) { if (c == '"') r += '"'; r += c; }
    return r + "\"";
}

std::string jsonString(const std::string& s) {
    std::string r = "\"";
    for (unsigned char c : s) {
        switch (c) {
            case '"': r += "\\\""; break;
            case '\\': r += "\\\\"; break;
            case '\n': r += "\\n"; break;
            case '\r': r += "\\r"; break;
            case '\t': r += "\\t"; break;
            default:
                if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); r += b; }
                else r += static_cast<char>(c);
        }
    }
    return r + "\"";
}

std::string num(double d) {
    char b[32];
    std::snprintf(b, sizeof b, "%.2f", d);
    return b;
}

}  // namespace

void writeFingerprintCsv(std::ostream& out, const std::vector<Fingerprint>& fps) {
    out << "ligand,residue,chain,resi,position,region";
    for (int c = 0; c < FpClassCount; c++) out << "," << fpClassName(c);
    out << "\n";
    for (const auto& fp : fps)
        for (const auto& e : fp.entries) {
            out << csvField(fp.ligandId) << "," << csvField(e.residue) << "," << csvField(e.chain) << "," << e.resi
                << (e.icode != ' ' ? std::string(1, e.icode) : "") << "," << csvField(e.position) << "," << csvField(e.region);
            for (int c = 0; c < FpClassCount; c++) out << "," << e.counts[c];
            out << "\n";
        }
}

void writeFingerprintJson(std::ostream& out, const std::vector<Fingerprint>& fps) {
    out << "[";
    for (size_t k = 0; k < fps.size(); k++) {
        const Fingerprint& fp = fps[k];
        out << (k ? ",\n " : "\n ") << "{\"ligand\": " << jsonString(fp.ligandId) << ", \"pocket\": " << jsonString(fp.pocketName)
            << ", \"classes\": [";
        for (int c = 0; c < FpClassCount; c++) out << (c ? ", " : "") << jsonString(fpClassName(c));
        out << "],\n  \"residues\": [";
        for (size_t i = 0; i < fp.entries.size(); i++) {
            const FingerprintEntry& e = fp.entries[i];
            out << (i ? ",\n   " : "\n   ") << "{\"residue\": " << jsonString(e.residue) << ", \"chain\": " << jsonString(e.chain)
                << ", \"resi\": " << e.resi << ", \"icode\": " << jsonString(e.icode == ' ' ? "" : std::string(1, e.icode))
                << ", \"position\": " << jsonString(e.position) << ", \"region\": " << jsonString(e.region)
                << ", \"in_pocket\": " << (e.inPocket ? "true" : "false") << ", \"counts\": {";
            for (int c = 0; c < FpClassCount; c++) out << (c ? ", " : "") << jsonString(fpClassName(c)) << ": " << e.counts[c];
            out << "}, \"min_distance\": {";
            bool first = true;
            for (int c = 0; c < FpClassCount; c++)
                if (e.counts[c] > 0) { out << (first ? "" : ", ") << jsonString(fpClassName(c)) << ": " << num(e.minDistance[c]); first = false; }
            out << "}}";
        }
        out << (fp.entries.empty() ? "" : "\n  ") << "]}";
    }
    out << (fps.empty() ? "" : "\n") << "]\n";
}

}  // namespace ligint
