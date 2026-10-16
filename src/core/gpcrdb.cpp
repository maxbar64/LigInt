// Pocket definitions for GPCRs from generic residue numbers: see include/ligint/gpcrdb.hpp.
#include "ligint/gpcrdb.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>

#include "internal.hpp"
#include "ligint/matrix.hpp"

namespace ligint {
namespace {

namespace fs = std::filesystem;

std::string trimmed(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::string lowered(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// The integer part of a generic number: the helix (1-7, 8 is helix 8) or the loop between two helices (12 = ICL1, 23 = ECL1, ...).
const std::map<int, const char*>& segmentNames() {
    static const std::map<int, const char*> names = {{1, "TM1"},  {2, "TM2"},  {3, "TM3"},  {4, "TM4"},  {5, "TM5"},  {6, "TM6"},  {7, "TM7"},
                                                     {8, "H8"},   {12, "ICL1"}, {23, "ECL1"}, {34, "ICL2"}, {45, "ECL2"}, {56, "ICL3"}, {67, "ECL3"}};
    return names;
}

// position along the receptor, to order a pocket: TM1, ICL1, TM2, ECL1, TM3, ICL2, TM4, ECL2, TM5, ICL3, TM6, ECL3, TM7, H8
int segmentRank(int ip) {
    static const std::map<int, int> rank = {{1, 0}, {12, 1}, {2, 2}, {23, 3}, {3, 4}, {34, 5}, {4, 6}, {45, 7}, {5, 8}, {56, 9}, {6, 10}, {67, 11}, {7, 12}, {8, 13}};
    const auto it = rank.find(ip);
    return it == rank.end() ? 99 : it->second;
}

const char* regionColor(const std::string& region) {
    static const std::map<std::string, const char*> colors = {
        {"TM1", "#ffa8a8"}, {"TM2", "#ffd8a8"}, {"TM3", "#ffec99"}, {"TM4", "#b2f2bb"}, {"TM5", "#99e9f2"}, {"TM6", "#a5d8ff"}, {"TM7", "#d0bfff"},
        {"H8", "#fcc2d7"},  {"ICL1", "#ced4da"}, {"ICL2", "#ced4da"}, {"ICL3", "#ced4da"}, {"ECL1", "#adb5bd"}, {"ECL2", "#adb5bd"}, {"ECL3", "#adb5bd"}};
    const auto it = colors.find(region);
    return it == colors.end() ? nullptr : it->second;
}

// the integer part and the digits after the separator of a canonical number
bool splitGeneric(const std::string& canonical, int& ip, std::string& digits) {
    std::size_t i = 0;
    while (i < canonical.size() && std::isdigit(static_cast<unsigned char>(canonical[i]))) ++i;
    if (i == 0 || i + 1 >= canonical.size()) return false;
    ip = std::atoi(canonical.substr(0, i).c_str());
    digits = canonical.substr(i + 1);
    return true;
}

// a canonical number comes before another in the order of the receptor
bool genericLess(const std::string& a, const std::string& b) {
    int ia = 0, ib = 0;
    std::string da, db;
    splitGeneric(a, ia, da);
    splitGeneric(b, ib, db);
    if (segmentRank(ia) != segmentRank(ib)) return segmentRank(ia) < segmentRank(ib);
    const std::size_t n = std::max(da.size(), db.size());
    da.resize(n, '0');
    db.resize(n, '0');
    return da < db;
}

std::vector<std::string> splitFields(const std::string& line, char sep) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (quoted) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
                else quoted = false;
            } else cur += c;
        } else if (c == '"' && trimmed(cur).empty()) {
            quoted = true;
            cur.clear();
        } else if (c == sep) {
            out.push_back(trimmed(cur));
            cur.clear();
        } else cur += c;
    }
    out.push_back(trimmed(cur));
    return out;
}

// A residue with the numbers that it has. gpcrdb and bw are canonical ("3x32", "3.32"), "" if the residue has none.
struct Numbered {
    std::string chain, res;
    int resi = 0;
    char icode = ' ';
    std::string gpcrdb, bw, segment;           // segment: only from a table (otherwise it comes from the number)
};

// ---------------------------------------------------------------------------------------------
// numbers in the B-factor column
// ---------------------------------------------------------------------------------------------
struct Res {
    std::string chain, res;
    int resi = 0;
    char icode = ' ';
    double ca = 0, n = 0;
    bool hasCa = false, hasN = false;
};

std::vector<Res> readResidues(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open the structure with generic numbers: " + path);
    while (in.peek() == '\n' || in.peek() == '\r' || in.peek() == ' ' || in.peek() == '\t') in.get();
    const int c = in.peek();
    const bool cif = (c == 'd' || c == 'D' || c == '#');                 // by content, as the analysis does
    std::vector<Res> residues;
    std::map<std::string, std::size_t> index;
    const detail::AtomSink sink = [&](detail::RawAtom& a) {
        if (a.het || (a.name != "CA" && a.name != "N")) return;
        const std::string key = a.chain + "|" + std::to_string(a.resi) + "|" + std::string(1, a.icode);
        auto it = index.find(key);
        if (it == index.end()) {
            Res r;
            r.chain = a.chain;
            r.res = a.res;
            r.resi = a.resi;
            r.icode = a.icode;
            index[key] = residues.size();
            residues.push_back(r);
            it = index.find(key);
        }
        Res& r = residues[it->second];
        if (a.name == "CA" && !r.hasCa) { r.ca = a.bf; r.hasCa = true; }          // the first alternate conformation
        if (a.name == "N" && !r.hasN) { r.n = a.bf; r.hasN = true; }
    };
    detail::CifStats stats;
    if (cif) detail::readMmcifAtoms(in, sink, stats);
    else detail::readPdbAtoms(in, sink);
    return residues;
}

struct ChainNumbers {
    std::string chain;
    std::vector<Numbered> items;               // the residues of the chain in file order
    std::size_t numbered = 0, numberedBw = 0;
    bool annotated = false, hasBw = false;
};

// The numbers of one chain, as GPCRdb writes them: the GPCRdb scheme in the B-factor of the CA atoms (the Ballesteros-Weinstein one in
// the N atoms). A residue without a number has 0.00. A bulge is the number of the preceding residue followed by "1": the raw output of
// GPCRdb marks it with a negative B-factor, and files that went through other programs repeat the number (a PDB B-factor has two decimals).
ChainNumbers numberChain(const std::string& chain, const std::vector<Res>& all, std::vector<std::string>& warnings) {
    ChainNumbers c;
    c.chain = chain;
    std::string prevBase, prevFinal;
    std::map<std::string, std::size_t> seen;
    std::vector<std::pair<int, int>> sequence, sequenceBw;           // (integer part, two decimals) of the residues that have a number, in order
    std::vector<bool> adjacent, adjacentBw;                          // whether the number follows a residue that has one
    bool lastNumbered = false, lastNumberedBw = false;
    for (const Res& r : all) {
        if (r.chain != chain) continue;
        Numbered n;
        n.chain = r.chain;
        n.res = r.res;
        n.resi = r.resi;
        n.icode = r.icode;
        bool has = false;
        if (r.hasCa && std::isfinite(r.ca) && r.ca != 0.0) {
            const long h = std::lround(std::fabs(r.ca) * 100.0);
            const int ip = static_cast<int>(h / 100), fp = static_cast<int>(h % 100);
            if (!(r.ca > 0 && h == 100) && segmentNames().count(ip)) {           // 1.00 is the filler of the programs that write B-factors
                char buf[24];
                std::snprintf(buf, sizeof buf, "%dx%02d", ip, fp);
                const std::string base = buf;
                std::string label;
                if (base == prevBase) label = prevFinal + "1";
                else if (r.ca < 0) label = base + "1";
                else label = base;
                if (seen.count(label)) {
                    warnings.push_back("the number " + label + " is on more than one residue of chain " + chain + ": the first is used");
                } else {
                    seen[label] = c.items.size();
                    n.gpcrdb = label;
                    has = true;
                    ++c.numbered;
                    sequence.push_back({ip, fp});
                    adjacent.push_back(lastNumbered);
                }
                prevBase = base;
                prevFinal = label;
            }
        }
        if (!has) { prevBase.clear(); prevFinal.clear(); }
        lastNumbered = has;
        bool hasBw = false;
        if (r.hasN && std::isfinite(r.n) && r.n > 0) {
            const long h = std::lround(r.n * 100.0);
            const int ip = static_cast<int>(h / 100), fp = static_cast<int>(h % 100);
            if (h != 100 && ip >= 1 && ip <= 8) {
                char buf[24];
                std::snprintf(buf, sizeof buf, "%d.%02d", ip, fp);
                n.bw = buf;
                hasBw = true;
                ++c.numberedBw;
                sequenceBw.push_back({ip, fp});
                adjacentBw.push_back(lastNumberedBw);
            }
        }
        lastNumberedBw = hasBw;
        c.items.push_back(n);
    }
    // Is this a numbering or the B-factors of a structure? Numbers grow along a helix; B-factors do not.
    auto nonDecreasing = [](const std::vector<std::pair<int, int>>& seq, const std::vector<bool>& adj, int& pairs) {
        int ok = 0;
        pairs = 0;
        for (std::size_t i = 1; i < seq.size(); ++i) {
            if (!adj[i] || seq[i].first != seq[i - 1].first) continue;
            ++pairs;
            ok += seq[i].second >= seq[i - 1].second;
        }
        return ok;
    };
    int pairs = 0;
    const int ok = nonDecreasing(sequence, adjacent, pairs);
    c.annotated = c.numbered >= 20 && c.numbered * 10 >= c.items.size() * 3 && (pairs < 8 || ok * 10 >= pairs * 9);
    int pairsBw = 0;
    const int okBw = nonDecreasing(sequenceBw, adjacentBw, pairsBw);
    c.hasBw = c.numberedBw >= 10 && (pairsBw < 8 || okBw * 10 >= pairsBw * 9);
    if (!c.hasBw)
        for (auto& n : c.items) n.bw.clear();
    return c;
}

// ---------------------------------------------------------------------------------------------
// a table
// ---------------------------------------------------------------------------------------------
std::string aa3(const std::string& code) {
    static const std::map<char, const char*> m = {{'A', "ALA"}, {'C', "CYS"}, {'D', "ASP"}, {'E', "GLU"}, {'F', "PHE"}, {'G', "GLY"}, {'H', "HIS"}, {'I', "ILE"}, {'K', "LYS"}, {'L', "LEU"},
                                                  {'M', "MET"}, {'N', "ASN"}, {'P', "PRO"}, {'Q', "GLN"}, {'R', "ARG"}, {'S', "SER"}, {'T', "THR"}, {'V', "VAL"}, {'W', "TRP"}, {'Y', "TYR"}};
    if (code.size() == 1) {
        const auto it = m.find(static_cast<char>(std::toupper(static_cast<unsigned char>(code[0]))));
        return it == m.end() ? std::string() : it->second;
    }
    std::string up = code;
    for (char& c : up) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return up.size() <= 3 ? up : std::string();
}

std::string tableColumn(const std::string& raw) {
    static const std::map<std::string, std::string> names = {
        {"chain", "chain"}, {"chain_id", "chain"}, {"auth_asym_id", "chain"},
        {"resi", "resi"}, {"resseq", "resi"}, {"residue", "resi"}, {"resnum", "resi"}, {"residue_number", "resi"}, {"sequence_number", "resi"}, {"number", "resi"}, {"pdb_resi", "resi"},
        {"gpcrdb", "gpcrdb"}, {"generic", "gpcrdb"}, {"generic_number", "gpcrdb"}, {"display_generic_number", "gpcrdb"}, {"gpcrdb_number", "gpcrdb"},
        {"bw", "bw"}, {"ballesteros", "bw"}, {"ballesteros_weinstein", "bw"}, {"bw_number", "bw"},
        {"segment", "segment"}, {"protein_segment", "segment"}, {"region", "segment"}, {"helix", "segment"},
        {"aa", "aa"}, {"amino_acid", "aa"}, {"resname", "aa"}, {"residue_name", "aa"}, {"res", "aa"}};
    const auto it = names.find(lowered(trimmed(raw)));
    return it == names.end() ? std::string() : it->second;
}

// A generic number in a table: "3x32", "3.32", or both as GPCRdb writes them, "3.32x32". Fills gpcrdb and/or bw; false if the text is not one.
bool tableNumber(const std::string& text, std::string& gpcrdb, std::string& bw) {
    const std::string t = trimmed(text);
    const std::size_t dot = t.find('.'), x = t.find_first_of("xX");
    if (dot != std::string::npos && x != std::string::npos && dot < x) {                 // both numbers: 3.32x32
        std::string c1, c2;
        char s1 = 0, s2 = 0;
        const std::string gp = t.substr(0, dot) + "x" + t.substr(x + 1);
        if (parseGenericNumber(t.substr(0, x), c1, s1) && s1 == '.' && parseGenericNumber(gp, c2, s2) && s2 == 'x') {
            bw = c1;
            gpcrdb = c2;
            return true;
        }
        return false;
    }
    std::string canonical;
    char scheme = 0;
    if (!parseGenericNumber(t, canonical, scheme)) return false;
    (scheme == 'x' ? gpcrdb : bw) = canonical;
    return true;
}

std::vector<Numbered> numberFromTable(const std::string& path, const std::string& defaultChain, std::vector<std::string>& warnings) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open the table of generic numbers: " + path);
    std::map<std::string, std::size_t> column;
    std::vector<Numbered> items;
    std::string line;
    char sep = '\t';
    bool header = false;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (lineNo == 1 && line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF && static_cast<unsigned char>(line[1]) == 0xBB && static_cast<unsigned char>(line[2]) == 0xBF) line.erase(0, 3);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string t = trimmed(line);
        if (t.empty() || t[0] == '#') continue;
        if (!header) {
            sep = line.find('\t') != std::string::npos ? '\t' : ',';
            const auto names = splitFields(line, sep);
            for (std::size_t i = 0; i < names.size(); ++i) {
                const std::string c = tableColumn(names[i]);
                if (!c.empty() && !column.count(c)) column[c] = i;
            }
            if (!column.count("resi") || (!column.count("gpcrdb") && !column.count("bw")))
                throw std::runtime_error("The table of generic numbers needs a header with a residue number column (resi, resseq...) and a column with the numbers "
                                         "(gpcrdb, generic... or bw)");
            header = true;
            continue;
        }
        const auto f = splitFields(line, sep);
        auto cell = [&](const char* name) {
            const auto it = column.find(name);
            return it == column.end() || it->second >= f.size() ? std::string() : f[it->second];
        };
        Numbered n;
        const std::string rs = cell("resi");
        if (rs.empty() || !parseLigandResidueNumber(rs, n.resi, n.icode))
            throw std::runtime_error("Line " + std::to_string(lineNo) + ": invalid residue number '" + rs + "' (expected e.g. 113 or 113A)");
        if (n.icode == 0) n.icode = ' ';
        n.chain = cell("chain").empty() ? defaultChain : cell("chain");
        n.res = aa3(cell("aa"));
        n.segment = cell("segment");
        for (const char* col : {"gpcrdb", "bw"}) {
            const std::string text = cell(col);
            if (text.empty() || text == "-" || text == "0") continue;
            std::string g, b;
            if (!tableNumber(text, g, b)) {
                warnings.push_back("line " + std::to_string(lineNo) + ": '" + text + "' is not a generic residue number and is ignored");
                continue;
            }
            if (!g.empty() && n.gpcrdb.empty()) n.gpcrdb = g;
            if (!b.empty() && n.bw.empty()) n.bw = b;
        }
        items.push_back(std::move(n));
    }
    if (!header) throw std::runtime_error("The table of generic numbers is empty");
    return items;
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// public functions
// ---------------------------------------------------------------------------------------------
bool parseGenericNumber(const std::string& text, std::string& canonical, char& scheme) {
    const std::string t = trimmed(text);
    std::size_t i = 0;
    while (i < t.size() && std::isdigit(static_cast<unsigned char>(t[i]))) ++i;
    if (i == 0 || i > 2 || i + 1 >= t.size()) return false;
    const char sep = t[i];
    if (sep != 'x' && sep != 'X' && sep != '.') return false;
    const std::string digits = t.substr(i + 1);
    if (digits.size() < 2 || digits.size() > 4) return false;
    for (char c : digits)
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    const int ip = std::atoi(t.substr(0, i).c_str());
    const bool bw = sep == '.';
    if (!segmentNames().count(ip)) return false;
    if (bw && (ip > 8 || digits.size() != 2)) return false;              // Ballesteros-Weinstein numbers: the helices, two decimals
    canonical = std::to_string(ip) + (bw ? "." : "x") + digits;
    scheme = bw ? '.' : 'x';
    return true;
}

std::string genericSegment(const std::string& canonical) {
    int ip = 0;
    std::string digits;
    if (!splitGeneric(canonical, ip, digits)) return std::string();
    const auto it = segmentNames().find(ip);
    return it == segmentNames().end() ? std::string() : it->second;
}

std::vector<std::string> parseGpcrdbPositions(const std::string& listOrFile) {
    std::string text = listOrFile;
    std::error_code ec;
    if (!listOrFile.empty() && fs::is_regular_file(fs::path(listOrFile), ec)) {
        std::ifstream f(listOrFile, std::ios::binary);
        if (!f) throw std::runtime_error("Cannot open the list of positions: " + listOrFile);
        std::ostringstream ss;
        ss << f.rdbuf();
        text = ss.str();
    }
    std::vector<std::string> out;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t' || c == '\r') c = ' ';
        std::istringstream words(line);
        std::string w;
        while (words >> w) {
            std::string canonical;
            char scheme = 0;
            if (!parseGenericNumber(w, canonical, scheme))
                throw std::runtime_error("Invalid generic residue number '" + w + "' (expected the helix and the position, e.g. 3x32 in the GPCRdb scheme or 3.32 in "
                                         "Ballesteros-Weinstein)");
            if (std::find(out.begin(), out.end(), canonical) == out.end()) out.push_back(canonical);
        }
    }
    if (out.empty()) throw std::runtime_error("No generic residue numbers given for the pocket");
    return out;
}

GpcrdbPocket generateGpcrdbPocket(const GpcrdbOptions& o) {
    if (o.positions.empty()) throw std::runtime_error("No generic residue numbers given for the pocket");
    std::vector<std::pair<std::string, char>> wanted;
    for (const auto& p : o.positions) {
        std::string canonical;
        char scheme = 0;
        if (!parseGenericNumber(p, canonical, scheme)) throw std::runtime_error("Invalid generic residue number '" + p + "'");
        if (std::none_of(wanted.begin(), wanted.end(), [&](const std::pair<std::string, char>& w) { return w.first == canonical; })) wanted.push_back({canonical, scheme});
    }
    GpcrdbPocket out;
    std::vector<Numbered> items;
    bool hasBw = false;
    std::string chain = o.chain;
    if (!o.mappingFile.empty()) {
        items = numberFromTable(o.mappingFile, o.chain, out.warnings);
        if (!o.chain.empty()) {
            std::vector<Numbered> mine;
            for (auto& n : items)
                if (n.chain == o.chain) mine.push_back(n);
            items.swap(mine);
        }
        for (const auto& n : items) hasBw = hasBw || !n.bw.empty();
        bool any = false;
        for (const auto& n : items) any = any || !n.gpcrdb.empty() || !n.bw.empty();
        if (!any) throw std::runtime_error("The table has no generic numbers" + (o.chain.empty() ? std::string() : " for chain " + o.chain));
    } else {
        if (o.numberedStructure.empty()) throw std::runtime_error("A structure annotated with generic numbers, or a table of them, is needed");
        const std::vector<Res> residues = readResidues(o.numberedStructure);
        std::vector<std::string> chainNames;
        for (const auto& r : residues)
            if (std::find(chainNames.begin(), chainNames.end(), r.chain) == chainNames.end()) chainNames.push_back(r.chain);
        if (chainNames.empty()) throw std::runtime_error("The structure has no protein residues with CA atoms, where GPCRdb writes the generic numbers: " + o.numberedStructure);
        std::vector<ChainNumbers> chains;
        for (const auto& cn : chainNames) {
            std::vector<std::string> w;
            chains.push_back(numberChain(cn, residues, w));
        }
        auto describe = [&]() {
            std::string text;
            for (const auto& c : chains) text += (text.empty() ? "" : ", ") + (c.chain.empty() ? std::string("(no id)") : c.chain) + ": " + std::to_string(c.numbered) + " of " + std::to_string(c.items.size());
            return text;
        };
        const ChainNumbers* chosen = nullptr;
        if (!o.chain.empty()) {
            for (const auto& c : chains)
                if (c.chain == o.chain) chosen = &c;
            if (!chosen) throw std::runtime_error("Chain " + o.chain + " is not in the structure (chains: " + describe() + ")");
        } else {
            std::vector<const ChainNumbers*> numbered;
            for (const auto& c : chains)
                if (c.annotated) numbered.push_back(&c);
            if (numbered.size() == 1) chosen = numbered[0];
            else if (numbered.size() > 1) {
                std::string names;
                for (auto* c : numbered) names += (names.empty() ? "" : ", ") + c->chain;
                throw std::runtime_error("Several chains have generic numbers (" + names + "): choose one with --gpcrdb-chain");
            }
        }
        if (!chosen || !chosen->annotated)
            throw std::runtime_error("The B-factor column of the structure does not have GPCRdb generic numbers (residues with a number, by chain: " + describe() + "). Annotate the PDB file "
                                     "with the tool \"Generic residue numbering (PDB)\" of GPCRdb (https://gpcrdb.org/structure/generic_numbering_index) and give the file that it returns");
        // the warnings of the chosen chain only (the other chains are not used)
        std::vector<std::string> w;
        const ChainNumbers again = numberChain(chosen->chain, residues, w);
        for (auto& x : w) out.warnings.push_back(x);
        items = again.items;
        hasBw = again.hasBw;
        chain = chosen->chain;
    }

    std::map<std::string, const Numbered*> byGpcrdb, byBw;
    for (const auto& n : items) {
        if (!n.gpcrdb.empty() && !byGpcrdb.count(n.gpcrdb)) byGpcrdb[n.gpcrdb] = &n;
        if (!n.bw.empty() && !byBw.count(n.bw)) byBw[n.bw] = &n;
    }
    bool needsBw = false;
    for (const auto& w : wanted) needsBw = needsBw || w.second == '.';
    if (needsBw && !hasBw)
        throw std::runtime_error("The Ballesteros-Weinstein numbers (with a point, e.g. 3.32) are not in this structure: GPCRdb writes them in the B-factor of the N atoms, or give "
                                 "a table with a bw column. The GPCRdb numbers (with an x, e.g. 3x32) are in the CA atoms. They are not the same near helix bulges and constrictions: "
                                 "for example 5.42 and 5x42 can be different residues");
    struct Found { std::string label; const Numbered* n; };
    std::vector<Found> found;
    for (const auto& w : wanted) {
        const auto& map = w.second == 'x' ? byGpcrdb : byBw;
        const auto it = map.find(w.first);
        if (it == map.end()) out.missing.push_back(w.first);
        else found.push_back({w.first, it->second});
    }
    if (found.empty()) {
        std::string list;
        for (const auto& w : wanted) list += (list.empty() ? "" : ", ") + w.first;
        throw std::runtime_error("None of the requested positions (" + list + ") has a residue in this structure" + (chain.empty() ? std::string() : " (chain " + chain + ")"));
    }
    std::stable_sort(found.begin(), found.end(), [](const Found& a, const Found& b) { return genericLess(a.label, b.label); });
    out.pocket.name = o.name.empty() ? "GPCR pocket (generic residue numbers" + (chain.empty() ? std::string() : ", chain " + chain) + ")" : o.name;
    std::set<std::string> regions;
    for (const auto& f : found) {
        PocketPosition p;
        p.chain = f.n->chain;
        p.resi = f.n->resi;
        p.icode = f.n->icode;
        p.label = f.label;
        p.region = f.n->segment.empty() ? genericSegment(f.label) : f.n->segment;
        p.resName = f.n->res;
        regions.insert(p.region);
        out.pocket.positions.push_back(p);
    }
    for (const auto& r : regions)
        if (const char* color = regionColor(r)) out.pocket.colors.push_back({r, color});
    if (!out.missing.empty()) {
        std::string list;
        for (const auto& m : out.missing) list += (list.empty() ? "" : ", ") + m;
        out.warnings.push_back("positions that this structure does not have (a missing or not numbered residue): " + list);
    }
    return out;
}

void saveGpcrdbPocket(const PocketDefinition& pocket, const std::string& filename) {
    std::ofstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write the pocket definition: " + filename);
    std::string name = pocket.name;
    for (char& c : name)
        if (c == '\n' || c == '\r') c = ' ';
    f << "@name " << name << "\n";
    for (const auto& c : pocket.colors) f << "@color " << c.first << " " << c.second << "\n";
    for (const auto& p : pocket.positions)
        f << p.chain << '\t' << p.resi << (p.icode == ' ' || p.icode == 0 ? std::string() : std::string(1, p.icode)) << '\t' << p.label << '\t' << p.region << '\t' << p.resName << '\n';
    f.flush();
    if (!f) throw std::runtime_error("Cannot write the pocket definition: " + filename);
}

}  // namespace ligint
