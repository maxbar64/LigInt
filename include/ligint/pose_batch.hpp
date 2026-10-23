// Summary of all imported ligand poses. No ranking or energy model is implied.
#pragma once
#include "ligint/ligint.hpp"
#include <iomanip>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace ligand_pose_batch {
struct Row {
    size_t index = 0, atoms = 0, favorable = 0, hbonds = 0, hydrophobic = 0;
    size_t pi = 0, salt = 0, water = 0, clashes = 0, repulsions = 0, covalent = 0, residues = 0;
    std::string name, error;
    std::map<std::string, std::set<std::string>> residueTypes;
    std::set<std::string> contactResidues; // unique favorable protein residues per pose
};
inline std::string csvQuote(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) { if (c == '"') out += "\"\""; else out += c; }
    return out + "\"";
}
inline std::string jsonQuote(const std::string& s) {
    std::ostringstream o;
    o << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') o << '\\' << char(c);
        else if (c == '\n') o << "\\n";
        else if (c == '\r') o << "\\r";
        else if (c == '\t') o << "\\t";
        else if (c < 32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else o << char(c);
    }
    o << '"';
    return o.str();
}
inline Row summarize(const ligint::Options& base, size_t idx, const ligint::Result& result) {
    Row row;
    row.index = idx + 1;
    row.name = base.templates[idx].name;
    for (const auto& a : base.templates[idx].at)
        if (a.el != "H" && a.el != "D") ++row.atoms;
    try {
        if (result.ligands.empty()) throw std::runtime_error("No ligand result");
        std::set<std::string> residues;
        for (const auto& h : result.ligands.front().hits) {
            if (h.favorable) ++row.favorable;
            if (h.type == ligint::type::HydrogenBond) ++row.hbonds;
            if (h.type == ligint::type::Hydrophobic || h.type == ligint::type::BackboneApolar) ++row.hydrophobic;
            if (h.type == ligint::type::PiStackingParallel || h.type == ligint::type::PiStackingT || h.type == ligint::type::CationPi) ++row.pi;
            if (h.type == ligint::type::SaltBridge) ++row.salt;
            if (h.type == ligint::type::WaterBridge) ++row.water;
            if (h.type == ligint::type::StericClash) ++row.clashes;
            if (h.type == ligint::type::ElectrostaticRepulsion) ++row.repulsions;
            if (h.type == ligint::type::PossibleCovalent) ++row.covalent;
            if (!h.residue.empty()) {
                residues.insert(h.residue);
                if (h.favorable) {
                    row.contactResidues.insert(h.residue);
                    row.residueTypes[h.residue].insert(h.type);
                }
            }
        }
        row.residues = residues.size();
    } catch (const std::exception& e) { row.error = e.what(); }
    return row;
}
inline Row analyze(const std::string& receptor, const ligint::Options& base, size_t idx) {
    Row row;
    row.index = idx + 1;
    row.name = base.templates[idx].name;
    try {
        ligint::Options opt = base;
        opt.externalLigandIndex = static_cast<int>(idx);
        return summarize(base, idx, ligint::analyzeFile(receptor, opt));
    } catch (const std::exception& e) {
        row.error = e.what();
        return row;
    }
}
// Exact residue/type fingerprint: geometry is intentionally not part of the signature.
inline std::string interactionSignature(const Row& row) {
    if (!row.error.empty()) return {};
    std::ostringstream signature;
    for (const auto& pair : row.residueTypes) {
        signature << pair.first.size() << ':' << pair.first << '=';
        for (const auto& type : pair.second)
            signature << type.size() << ':' << type << ';';
        signature << '|';
    }
    return signature.str();
}
inline std::vector<size_t> exactFingerprintGroups(const std::vector<Row>& rows) {
    std::map<std::string, size_t> known;
    std::vector<size_t> groups;
    groups.reserve(rows.size());
    for (const Row& row : rows) {
        if (!row.error.empty()) { groups.push_back(0); continue; }
        const std::string signature = interactionSignature(row);
        auto inserted = known.emplace(signature, known.size() + 1);
        groups.push_back(inserted.first->second);
    }
    return groups;
}
inline void writeResidueMatrix(std::ostream& o, const std::vector<Row>& rows, const std::string& format) {
    std::map<std::string, size_t> counts;
    size_t successful = 0;
    for (const auto& row : rows) {
        if (!row.error.empty()) continue;
        ++successful;
        for (const auto& pair : row.residueTypes) ++counts[pair.first];
    }
    if (format == "csv") {
        o << "residue,pose_count,frequency_percent";
        for (const auto& row : rows) o << ',' << csvQuote("pose_" + std::to_string(row.index));
        o << '\n';
        for (const auto& [residue, count] : counts) {
            o << csvQuote(residue) << ',' << count << ',' << (successful ? 100.0 * count / successful : 0.0);
            for (const auto& row : rows) {
                std::string types;
                if (!row.error.empty()) types = "ERROR";
                else {
                    auto found = row.residueTypes.find(residue);
                    if (found != row.residueTypes.end())
                        for (const auto& type : found->second) {
                            if (!types.empty()) types += "; ";
                            types += type;
                        }
                }
                o << ',' << csvQuote(types);
            }
            o << '\n';
        }
    } else {
        o << "{\n  \"successful_poses\":" << successful << ",\n  \"residues\":[\n";
        bool firstResidue = true;
        for (const auto& [residue, count] : counts) {
            if (!firstResidue) o << ",\n";
            firstResidue = false;
            o << "    {\"residue\":" << jsonQuote(residue) << ",\"pose_count\":" << count
              << ",\"frequency_percent\":" << (successful ? 100.0 * count / successful : 0.0) << ",\"poses\":[";
            for (size_t i = 0; i < rows.size(); ++i) {
                if (i) o << ',';
                const auto& row = rows[i];
                o << "{\"index\":" << row.index << ",\"status\":" << jsonQuote(row.error.empty() ? "ok" : "error") << ",\"types\":[";
                auto found = row.residueTypes.find(residue);
                if (row.error.empty() && found != row.residueTypes.end()) {
                    bool firstType = true;
                    for (const auto& type : found->second) {
                        if (!firstType) o << ',';
                        firstType = false;
                        o << jsonQuote(type);
                    }
                }
                o << "]}";
            }
            o << "]}";
        }
        o << "\n  ]\n}\n";
    }
}
inline void write(std::ostream& o, const std::vector<Row>& rows, const std::string& format) {
    if (format == "csv") {
        o << "pose_index,pose_name,status,error,heavy_atoms,favorable_contacts,hbonds,hydrophobic,pi,salt_bridges,water_bridges,clashes,repulsions,possible_covalent,interacting_residues\n";
        for (const auto& r : rows)
            o << r.index << ',' << csvQuote(r.name) << ',' << (r.error.empty() ? "ok" : "error") << ',' << csvQuote(r.error)
              << ',' << r.atoms << ',' << r.favorable << ',' << r.hbonds << ',' << r.hydrophobic << ',' << r.pi << ',' << r.salt << ',' << r.water
              << ',' << r.clashes << ',' << r.repulsions << ',' << r.covalent << ',' << r.residues << '\n';
    } else {
        o << "{\n  \"poses\": [\n";
        for (size_t i = 0; i < rows.size(); ++i) {
            const auto& r = rows[i];
            if (i) o << ",\n";
            o << "    {\"pose_index\":" << r.index << ",\"pose_name\":" << jsonQuote(r.name) << ",\"status\":" << jsonQuote(r.error.empty() ? "ok" : "error")
              << ",\"error\":" << jsonQuote(r.error) << ",\"heavy_atoms\":" << r.atoms << ",\"favorable_contacts\":" << r.favorable
              << ",\"hbonds\":" << r.hbonds << ",\"hydrophobic\":" << r.hydrophobic << ",\"pi\":" << r.pi << ",\"salt_bridges\":" << r.salt
              << ",\"water_bridges\":" << r.water << ",\"clashes\":" << r.clashes << ",\"repulsions\":" << r.repulsions
              << ",\"possible_covalent\":" << r.covalent << ",\"interacting_residues\":" << r.residues << '}';
        }
        o << "\n  ]\n}\n";
    }
}
} // namespace ligand_pose_batch
