// Interaction matrix: see include/ligint/matrix.hpp.
#include "ligint/matrix.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

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

// Splits a line of a TSV/CSV file; a field in double quotes may contain the separator ("" is a quote).
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

std::string canonicalColumn(const std::string& raw) {
    static const std::map<std::string, std::string> names = {
        {"structure", "structure"}, {"file", "structure"}, {"path", "structure"}, {"structure_file", "structure"}, {"pdb", "structure"},
        {"ligand", "ligand"}, {"lig", "ligand"}, {"resname", "ligand"}, {"ligand_name", "ligand"},
        {"chain", "chain"}, {"ligand_chain", "chain"},
        {"resseq", "resseq"}, {"resi", "resseq"}, {"residue", "resseq"}, {"resnum", "resseq"}, {"number", "resseq"}, {"ligand_resi", "resseq"},
        {"label", "label"}, {"name", "label"}, {"id", "label"}, {"entry", "label"},
        {"pocket", "pocket"}, {"pocket_file", "pocket"}, {"definition", "pocket"},
        {"altloc", "altloc"}, {"alt", "altloc"}, {"conformation", "altloc"}};
    const auto it = names.find(lowered(trimmed(raw)));
    return it == names.end() ? std::string() : it->second;
}

std::string resolvePath(const std::string& p, const std::string& baseDir) {
    if (p.empty() || baseDir.empty()) return p;
    const fs::path q(p);
    if (q.is_absolute()) return p;
    return (fs::path(baseDir) / q).string();
}

std::string csvEscape(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += '"';
        out += c;
    }
    return out + "\"";
}

std::string jsonEscape(const std::string& s) {
    std::string out;
    for (unsigned char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof buf, "\\u%04x", c);
                    out += buf;
                } else out += static_cast<char>(c);
        }
    }
    return out;
}

std::string xmlEscape(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c;
        }
    }
    return out;
}

// the first n bytes, without cutting a UTF-8 character
std::string prefixOf(const std::string& s, std::size_t n) {
    if (s.size() <= n) return s;
    std::size_t k = n;
    while (k > 0 && (static_cast<unsigned char>(s[k]) & 0xC0) == 0x80) --k;
    return s.substr(0, k);
}

std::string truncated(const std::string& s, std::size_t n) { return s.size() <= n ? s : prefixOf(s, n) + "..."; }

std::string hetList(const Result& r) {
    std::string out;
    for (const auto& kv : r.hetCounts) out += (out.empty() ? "" : " ") + kv.first + "(" + std::to_string(kv.second) + ")";
    return out.empty() ? std::string("no ligands") : out;
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// the list
// ---------------------------------------------------------------------------------------------
bool parseLigandResidueNumber(const std::string& rs, int& resi, char& icode) {
    std::size_t k = 0;
    if (k < rs.size() && (rs[k] == '-' || rs[k] == '+')) ++k;
    const std::size_t d0 = k;
    while (k < rs.size() && std::isdigit(static_cast<unsigned char>(rs[k]))) ++k;
    const std::string suffix = rs.substr(k);
    const bool okSuffix = suffix.empty() || suffix == "-" || (suffix.size() == 1 && std::isalpha(static_cast<unsigned char>(suffix[0])));
    if (k == d0 || k - d0 > 6 || !okSuffix) return false;
    resi = std::atoi(rs.substr(0, k).c_str());
    icode = suffix.empty() ? 0 : (suffix == "-" ? ' ' : suffix[0]);
    return true;
}

std::vector<MatrixEntry> parseMatrixList(std::istream& in, const std::string& baseDir) {
    std::vector<MatrixEntry> entries;
    std::map<std::string, std::size_t> column;      // canonical name -> index in the line
    std::string line;
    char sep = '\t';
    bool header = false;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (lineNo == 1 && line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF && static_cast<unsigned char>(line[1]) == 0xBB &&
            static_cast<unsigned char>(line[2]) == 0xBF)
            line.erase(0, 3);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string t = trimmed(line);
        if (t.empty() || t[0] == '#') continue;
        if (!header) {
            sep = line.find('\t') != std::string::npos ? '\t' : ',';
            const auto names = splitFields(line, sep);
            for (std::size_t i = 0; i < names.size(); ++i) {
                const std::string c = canonicalColumn(names[i]);
                if (!c.empty() && !column.count(c)) column[c] = i;
            }
            if (!column.count("structure"))
                throw std::runtime_error("Line " + std::to_string(lineNo) + ": the list needs a header line with a 'structure' column "
                                         "(columns: structure, ligand, chain, resseq, label, pocket, altloc)");
            header = true;
            continue;
        }
        const auto f = splitFields(line, sep);
        auto cell = [&](const char* name) {
            const auto it = column.find(name);
            return it == column.end() || it->second >= f.size() ? std::string() : f[it->second];
        };
        const std::string where = "Line " + std::to_string(lineNo) + ": ";
        MatrixEntry e;
        e.line = lineNo;
        e.structure = cell("structure");
        if (e.structure.empty()) throw std::runtime_error(where + "the structure is empty");
        e.structure = resolvePath(e.structure, baseDir);
        e.ligand = cell("ligand");
        e.chain = cell("chain");
        const std::string rs = cell("resseq");
        if (!rs.empty() && !parseLigandResidueNumber(rs, e.resi, e.icode))
            throw std::runtime_error(where + "invalid residue number '" + rs + "' (expected e.g. 201, 201A or 201-)");
        e.label = cell("label");
        e.pocket = resolvePath(cell("pocket"), baseDir);
        const std::string alt = cell("altloc");
        if (!alt.empty()) {
            if (alt.size() != 1 || alt[0] == ' ') throw std::runtime_error(where + "altloc expects a single character, e.g. A or B");
            e.altloc = alt[0];
        }
        entries.push_back(std::move(e));
        if (entries.size() > 100000) throw std::runtime_error("The list has more than 100000 entries");
    }
    if (!header) throw std::runtime_error("The list is empty (it needs a header line and one line for each ligand)");
    if (entries.empty()) throw std::runtime_error("The list has a header but no entries");
    return entries;
}

std::vector<MatrixEntry> readMatrixList(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open the matrix list: " + path);
    return parseMatrixList(in, fs::path(path).parent_path().string());
}

// ---------------------------------------------------------------------------------------------
// classes
// ---------------------------------------------------------------------------------------------
std::vector<int> defaultMatrixClasses() { return {FpHydrophobic, FpFaceToFace, FpEdgeToFace, FpHBondDonor, FpHBondAcceptor, FpIonicPositive, FpIonicNegative}; }

bool parseMatrixClasses(const std::string& text, std::vector<int>& out, std::string& error) {
    out.clear();
    const std::string t = lowered(trimmed(text));
    if (t == "klifs") { out = defaultMatrixClasses(); return true; }
    if (t == "all") {
        for (int c = 0; c < FpClassCount; ++c) out.push_back(c);
        return true;
    }
    std::set<int> chosen;
    std::stringstream ss(t);
    std::string name;
    while (std::getline(ss, name, ',')) {
        name = trimmed(name);
        if (name.empty()) continue;
        int found = -1;
        for (int c = 0; c < FpClassCount; ++c)
            if (name == fpClassName(c)) found = c;
        if (found < 0) {
            error = "unknown interaction class '" + name + "' (use klifs, all or names such as hydrophobic,hbond_donor)";
            return false;
        }
        chosen.insert(found);
    }
    if (chosen.empty()) {
        error = "no interaction class given (use klifs, all or a comma-separated list of names)";
        return false;
    }
    out.assign(chosen.begin(), chosen.end());
    return true;
}

int MatrixRow::bitCount() const {
    int n = 0;
    for (unsigned char b : bits) n += b != 0;
    return n;
}

// ---------------------------------------------------------------------------------------------
// building
// ---------------------------------------------------------------------------------------------
InteractionMatrix makeMatrix(const std::vector<MatrixEntry>& entries, const MatrixOptions& options) {
    if (entries.empty()) throw std::runtime_error("The list has no entries");
    InteractionMatrix m;
    m.classes = options.classes.empty() ? defaultMatrixClasses() : options.classes;
    for (int c : m.classes)
        if (c < 0 || c >= FpClassCount) throw std::runtime_error("Invalid interaction class");
    const std::size_t nc = m.classes.size();

    // all the entries use a pocket definition, or none does: the columns must mean the same thing for every row
    std::size_t withPocket = 0;
    for (const auto& e : entries) withPocket += !(e.pocket.empty() && options.pocket.empty());
    if (withPocket != 0 && withPocket != entries.size()) {
        for (std::size_t i = 0; i < entries.size(); ++i)
            if (entries[i].pocket.empty() && options.pocket.empty())
                throw std::runtime_error("Entry " + std::to_string(i + 1) + (entries[i].line > 0 ? " (line " + std::to_string(entries[i].line) + ")" : std::string()) +
                                         " has no pocket definition while others have one: give every entry a pocket (or a common one), or none");
    }
    m.pocketMode = withPocket != 0;

    struct Tmp {
        std::string label, structure, ligandId;
        std::vector<FingerprintEntry> items;
    };
    std::vector<Tmp> tmp;
    std::set<std::string> labels;
    std::map<std::string, PocketDefinition> pockets;

    for (std::size_t idx = 0; idx < entries.size(); ++idx) {
        const MatrixEntry& e = entries[idx];
        const std::string name = e.label.empty() ? fs::path(e.structure).filename().string() : e.label;
        const std::string where = "Entry " + std::to_string(idx + 1) + (e.line > 0 ? " (line " + std::to_string(e.line) + ", " : " (") + name + "): ";
        Options o = options.analysis;
        o.ligName = e.ligand;
        o.chain = e.chain;
        o.resi = e.resi;
        o.icode = e.icode;
        if (e.altloc) o.altloc = e.altloc;
        Result r;
        try {
            r = analyzeFile(e.structure, o);
            if (e.ligand.empty()) throw std::runtime_error("the ligand is required; the structure has: " + hetList(r));
        } catch (const std::runtime_error& ex) {
            throw std::runtime_error(where + ex.what());
        }
        if (r.listingOnly) throw std::runtime_error(where + "ligand '" + e.ligand + "' not found; the structure has: " + hetList(r));
        if (r.ligands.size() != 1) {
            std::string ids;
            for (const auto& l : r.ligands) ids += (ids.empty() ? "" : ", ") + l.id;
            throw std::runtime_error(where + "'" + e.ligand + "' matches " + std::to_string(r.ligands.size()) + " ligands (" + ids +
                                     "): add the chain and the residue number to choose one");
        }
        const LigandResult& lig = r.ligands[0];
        const std::string label = e.label.empty() ? fs::path(e.structure).stem().string() + " " + lig.id : e.label;
        if (!labels.insert(label).second) throw std::runtime_error(where + "the label '" + label + "' is used by another entry: give each row its own label");
        for (const auto& w : r.warnings) m.warnings.push_back(label + ": " + w);

        const PocketDefinition* pocket = nullptr;
        if (m.pocketMode) {
            const std::string& path = e.pocket.empty() ? options.pocket : e.pocket;
            auto it = pockets.find(path);
            if (it == pockets.end()) {
                PocketDefinition pd;
                std::vector<std::string> pw;
                if (!loadPocketDefinition(path, pd, pw)) throw std::runtime_error(where + "cannot open the pocket definition " + path);
                for (const auto& w : pw) m.warnings.push_back(path + ": " + w);
                it = pockets.emplace(path, std::move(pd)).first;
            }
            pocket = &it->second;
        }
        Fingerprint fp = makeFingerprint(lig, pocket);
        for (const auto& w : fp.warnings) m.warnings.push_back(label + ": " + w);
        for (const auto& rc : fp.regionColors) {
            bool known = false;
            for (const auto& have : m.regionColors) known = known || have.first == rc.first;
            if (!known) m.regionColors.push_back(rc);
        }
        Tmp t;
        t.label = label;
        t.structure = e.structure;
        t.ligandId = lig.id;
        t.items = std::move(fp.entries);
        tmp.push_back(std::move(t));
    }

    // the columns
    std::map<std::string, std::size_t> index;
    if (m.pocketMode) {
        // The positions of the pockets, in order. A position that only some pockets have goes right after the one that precedes it
        // in its pocket (so the numeric order of KLIFS positions survives pockets with missing residues).
        std::vector<std::string> order;
        std::map<std::string, std::string> regionOf;
        for (const auto& t : tmp) {
            std::vector<std::string> mine;
            std::vector<std::string> regions;
            for (const auto& it : t.items)
                if (it.inPocket && !it.position.empty() && std::find(mine.begin(), mine.end(), it.position) == mine.end()) {
                    mine.push_back(it.position);
                    regions.push_back(it.region);
                }
            for (std::size_t i = 0; i < mine.size(); ++i) {
                const auto known = regionOf.find(mine[i]);
                if (known != regionOf.end()) {
                    if (!regions[i].empty() && known->second != regions[i])
                        m.warnings.push_back("position " + mine[i] + " is in the region " + known->second + " in one pocket definition and in " + regions[i] +
                                             " in another: the first is used");
                    continue;
                }
                std::size_t at = 0;
                if (i > 0) at = static_cast<std::size_t>(std::find(order.begin(), order.end(), mine[i - 1]) - order.begin()) + 1;
                order.insert(order.begin() + static_cast<std::ptrdiff_t>(at), mine[i]);
                regionOf[mine[i]] = regions[i];
            }
        }
        for (const auto& label : order) {
            index[label] = m.positions.size();
            m.positions.push_back({label, regionOf[label]});
        }
    } else {
        std::map<std::tuple<std::string, int, char, std::string>, std::string> seen;
        for (const auto& t : tmp)
            for (const auto& it : t.items)
                if (it.any()) seen[std::make_tuple(it.chain, it.resi, it.icode, it.resName)] = it.residue;
        for (const auto& kv : seen) {
            if (index.count(kv.second)) continue;
            index[kv.second] = m.positions.size();
            m.positions.push_back({kv.second, std::string()});
        }
    }

    // the rows
    const std::size_t np = m.positions.size();
    for (const auto& t : tmp) {
        MatrixRow row;
        row.label = t.label;
        row.structure = t.structure;
        row.ligandId = t.ligandId;
        row.residues.assign(np, std::string());
        row.bits.assign(np * nc, 0);
        std::vector<std::string> outside;
        for (const auto& it : t.items) {
            bool selected = false;
            for (int c : m.classes) selected = selected || it.counts[static_cast<std::size_t>(c)] > 0;
            std::string key;
            if (m.pocketMode) {
                if (!it.inPocket || it.position.empty()) {
                    if (selected) outside.push_back(it.residue);
                    continue;
                }
                key = it.position;
            } else {
                if (!it.any()) continue;
                key = it.residue;
            }
            const auto at = index.find(key);
            if (at == index.end()) continue;
            row.residues[at->second] = it.residue;
            for (std::size_t ci = 0; ci < nc; ++ci)
                if (it.counts[static_cast<std::size_t>(m.classes[ci])] > 0) row.bits[at->second * nc + ci] = 1;
        }
        if (!outside.empty()) {
            std::string list;
            for (std::size_t i = 0; i < outside.size() && i < 5; ++i) list += (i ? ", " : "") + outside[i];
            m.warnings.push_back(t.label + ": " + std::to_string(outside.size()) + " residue(s) in contact are not in the pocket definition and are not in the matrix: " + list +
                                 (outside.size() > 5 ? ", ..." : ""));
        }
        m.rows.push_back(std::move(row));
    }

    if (options.hideEmptyPositions) {
        std::vector<std::size_t> keep;
        for (std::size_t p = 0; p < np; ++p) {
            bool any = false;
            for (const auto& row : m.rows)
                for (std::size_t ci = 0; ci < nc; ++ci) any = any || row.bits[p * nc + ci] != 0;
            if (any) keep.push_back(p);
        }
        std::vector<MatrixPosition> positions;
        for (std::size_t p : keep) positions.push_back(m.positions[p]);
        for (auto& row : m.rows) {
            std::vector<unsigned char> bits;
            std::vector<std::string> residues;
            for (std::size_t p : keep) {
                residues.push_back(row.residues[p]);
                for (std::size_t ci = 0; ci < nc; ++ci) bits.push_back(row.bits[p * nc + ci]);
            }
            row.bits = std::move(bits);
            row.residues = std::move(residues);
        }
        m.positions = std::move(positions);
    }
    return m;
}

// ---------------------------------------------------------------------------------------------
// CSV and JSON
// ---------------------------------------------------------------------------------------------
void writeMatrixCsv(std::ostream& out, const InteractionMatrix& m) {
    out << "label,structure,ligand";
    for (const auto& p : m.positions)
        for (int c : m.classes) out << ',' << csvEscape(p.label + ":" + fpClassName(c));
    out << '\n';
    for (const auto& row : m.rows) {
        out << csvEscape(row.label) << ',' << csvEscape(row.structure) << ',' << csvEscape(row.ligandId);
        for (unsigned char b : row.bits) out << ',' << (b ? '1' : '0');
        out << '\n';
    }
}

void writeMatrixJson(std::ostream& out, const InteractionMatrix& m) {
    out << "{\n  \"classes\": [";
    for (std::size_t i = 0; i < m.classes.size(); ++i) out << (i ? ", " : "") << '"' << fpClassName(m.classes[i]) << '"';
    out << "],\n  \"pocket\": " << (m.pocketMode ? "true" : "false") << ",\n  \"positions\": [";
    for (std::size_t i = 0; i < m.positions.size(); ++i)
        out << (i ? ", " : "") << "{\"position\": \"" << jsonEscape(m.positions[i].label) << "\", \"region\": \"" << jsonEscape(m.positions[i].region) << "\"}";
    out << "],\n  \"rows\": [";
    for (std::size_t r = 0; r < m.rows.size(); ++r) {
        const MatrixRow& row = m.rows[r];
        out << (r ? "," : "") << "\n    {\"label\": \"" << jsonEscape(row.label) << "\", \"structure\": \"" << jsonEscape(row.structure) << "\", \"ligand\": \""
            << jsonEscape(row.ligandId) << "\",\n     \"residues\": [";
        for (std::size_t i = 0; i < row.residues.size(); ++i) out << (i ? ", " : "") << '"' << jsonEscape(row.residues[i]) << '"';
        out << "],\n     \"bits\": [";
        for (std::size_t i = 0; i < row.bits.size(); ++i) out << (i ? "," : "") << (row.bits[i] ? '1' : '0');
        out << "]}";
    }
    out << (m.rows.empty() ? "" : "\n  ") << "],\n  \"warnings\": [";
    for (std::size_t i = 0; i < m.warnings.size(); ++i) out << (i ? ", " : "") << '"' << jsonEscape(m.warnings[i]) << '"';
    out << "]\n}\n";
}

// ---------------------------------------------------------------------------------------------
// SVG heat map
// ---------------------------------------------------------------------------------------------
std::string matrixSvg(const InteractionMatrix& m, const MatrixSvgOptions& opt) {
    // one distinct color for each class (distinguishable in the neighbouring stripes of a block)
    static const char* const color[FpClassCount] = {"#8c8c8c", "#cc79a7", "#e69f00", "#0072b2", "#d55e00", "#56b4e9", "#009e73",
                                                    "#6a3d9a", "#b15928", "#222222", "#a6cee3"};
    static const char* const bandColor[] = {"#a5d8ff", "#ffd8a8", "#d0bfff", "#b2f2bb", "#ffc9c9", "#fff3bf", "#99e9f2", "#fcc2d7"};
    const int stripe = 3, stripeH = 12, rowPitch = 15, gap = 3;
    const int nc = static_cast<int>(m.classes.size());
    const int np = static_cast<int>(m.positions.size());
    const int blockW = nc * stripe + gap;

    std::size_t longestLabel = 0;
    for (const auto& r : m.rows) longestLabel = std::max(longestLabel, std::min<std::size_t>(r.label.size(), 44));
    const int leftW = std::max(90, static_cast<int>(longestLabel) * 6 + 18);
    std::size_t longestPosition = 1;
    for (const auto& p : m.positions) longestPosition = std::max(longestPosition, p.label.size());
    const int posH = std::min(80, static_cast<int>(longestPosition) * 5 + 8);
    bool hasRegions = false;
    for (const auto& p : m.positions) hasRegions = hasRegions || !p.region.empty();
    const int titleH = opt.title.empty() ? 0 : 24;
    const int bandH = hasRegions ? 16 : 0;
    const int y0 = titleH + bandH + posH + 4;
    const int rows = static_cast<int>(m.rows.size());
    const int width = leftW + np * blockW + 24;
    // the legend wraps to the width of the figure: lay it out first to know its height
    struct LegendItem { int x, line, cls; };
    std::vector<LegendItem> legend;
    int legendLines = 0;
    if (opt.showLegend) {
        int x = leftW, line = 0;
        for (int cls : m.classes) {
            const int w = 16 + static_cast<int>(std::string(fpClassTitle(cls)).size()) * 6 + 14;
            if (x + w > width && x > leftW) { x = leftW; ++line; }
            legend.push_back({x, line, cls});
            x += w;
        }
        legendLines = line + 1;
    }
    const int height = y0 + rows * rowPitch + 12 + (opt.showLegend ? legendLines * 14 + 12 : 0);

    std::ostringstream s;
    s << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << ' ' << height
      << "\" font-family=\"sans-serif\">\n<rect width=\"" << width << "\" height=\"" << height << "\" fill=\"#ffffff\"/>\n";
    if (!opt.title.empty()) s << "<text x=\"" << leftW << "\" y=\"16\" font-size=\"13\" font-weight=\"bold\" fill=\"#212529\">" << xmlEscape(opt.title) << "</text>\n";

    // the shading of the columns, the regions and the position labels
    for (int p = 0; p < np; ++p)
        if (p % 2 == 0) s << "<rect x=\"" << leftW + p * blockW << "\" y=\"" << y0 << "\" width=\"" << blockW << "\" height=\"" << rows * rowPitch << "\" fill=\"#f1f3f5\"/>\n";
    if (hasRegions) {
        int start = 0, regionIndex = 0;
        std::map<std::string, int> seen;
        while (start < np) {
            int end = start;
            while (end + 1 < np && m.positions[static_cast<std::size_t>(end + 1)].region == m.positions[static_cast<std::size_t>(start)].region) ++end;
            const std::string& region = m.positions[static_cast<std::size_t>(start)].region;
            std::string fill;
            for (const auto& rc : m.regionColors)
                if (rc.first == region) fill = rc.second;
            if (fill.empty()) {
                if (!seen.count(region)) seen[region] = regionIndex++;
                fill = bandColor[seen[region] % 8];
            }
            const int x = leftW + start * blockW, w = (end - start + 1) * blockW - 1;
            s << "<rect x=\"" << x << "\" y=\"" << titleH << "\" width=\"" << w << "\" height=\"" << bandH - 2 << "\" fill=\"" << xmlEscape(fill) << "\"/>\n";
            if (!region.empty()) {
                // the name is cut to the width of its band (the full name is in the tooltip)
                const std::size_t fits = static_cast<std::size_t>(std::max(1, (w - 2) / 5));
                // the tooltip is a <title> beside the text: inside a <text> the Qt renderer would draw it as a second text
                s << "<g><title>" << xmlEscape(region) << "</title><text x=\"" << x + w / 2 << "\" y=\"" << titleH + bandH - 5
                  << "\" font-size=\"8\" text-anchor=\"middle\" fill=\"#212529\">" << xmlEscape(prefixOf(region, fits)) << "</text></g>\n";
            }
            start = end + 1;
        }
    }
    for (int p = 0; p < np; ++p)
        s << "<text transform=\"translate(" << leftW + p * blockW + blockW / 2 + 2 << ' ' << y0 - 3 << ") rotate(-90)\" font-size=\"8\" fill=\"#495057\">"
          << xmlEscape(m.positions[static_cast<std::size_t>(p)].label) << "</text>\n";

    // the rows
    for (int r = 0; r < rows; ++r) {
        const MatrixRow& row = m.rows[static_cast<std::size_t>(r)];
        const int y = y0 + r * rowPitch;
        s << "<line x1=\"" << leftW << "\" y1=\"" << y + rowPitch << "\" x2=\"" << leftW + np * blockW << "\" y2=\"" << y + rowPitch << "\" stroke=\"#dee2e6\" stroke-width=\"0.5\"/>\n";
        s << "<text x=\"" << leftW - 6 << "\" y=\"" << y + 10 << "\" font-size=\"9\" text-anchor=\"end\" fill=\"#212529\">" << xmlEscape(truncated(row.label, 44)) << "</text>\n";
        for (int p = 0; p < np; ++p)
            for (int c = 0; c < nc; ++c) {
                const std::size_t at = static_cast<std::size_t>(p * nc + c);
                if (!row.bits[at]) continue;
                const int cls = m.classes[static_cast<std::size_t>(c)];
                s << "<rect x=\"" << leftW + p * blockW + c * stripe << "\" y=\"" << y + 1 << "\" width=\"" << stripe << "\" height=\"" << stripeH << "\" fill=\"" << color[cls]
                  << "\"><title>" << xmlEscape(row.label + ": " + m.positions[static_cast<std::size_t>(p)].label + " " + fpClassTitle(cls) +
                                              (row.residues[static_cast<std::size_t>(p)].empty() ? std::string() : " (" + row.residues[static_cast<std::size_t>(p)] + ")"))
                  << "</title></rect>\n";
            }
    }

    for (const auto& item : legend) {
        const int y = y0 + rows * rowPitch + 26 + item.line * 14;
        s << "<rect x=\"" << item.x << "\" y=\"" << y - 9 << "\" width=\"10\" height=\"10\" fill=\"" << color[item.cls] << "\"/><text x=\"" << item.x + 14 << "\" y=\"" << y
          << "\" font-size=\"9\" fill=\"#212529\">" << xmlEscape(fpClassTitle(item.cls)) << "</text>\n";
    }
    s << "</svg>\n";
    return s.str();
}

}  // namespace ligint
