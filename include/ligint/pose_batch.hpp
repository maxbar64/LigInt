// Summary of all imported ligand poses. No ranking or energy model is implied.
#pragma once
#include "ligint/ligint.hpp"
#include <iomanip>
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
inline Row analyze(const std::string& receptor, const ligint::Options& base, size_t idx) {
    Row row;
    row.index = idx + 1;
    row.name = base.templates[idx].name;
    for (const auto& a : base.templates[idx].at)
        if (a.el != "H" && a.el != "D") ++row.atoms;
    try {
        ligint::Options opt = base;
        opt.externalLigandIndex = static_cast<int>(idx);
        const auto result = ligint::analyzeFile(receptor, opt);
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
                if (h.favorable) row.contactResidues.insert(h.residue);
            }
        }
        row.residues = residues.size();
    } catch (const std::exception& e) { row.error = e.what(); }
    return row;
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
