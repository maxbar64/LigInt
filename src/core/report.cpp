// Text reports: readable table and CSV.
#include "ligint/ligint.hpp"

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace ligint {
namespace {

std::string fmt(const char* f, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

void printListing(std::ostream& out, const Result& r) {
    out << "Ligands (HETATM) present:";
    for (auto& kv : r.hetCounts) out << " " << kv.first << "(" << kv.second << ")";
    out << "\n";
    if (!r.notFound.empty()) out << "No ligand '" << r.notFound << "' found.\n";
}

}  // namespace

void printText(std::ostream& out, const Result& r) {
    if (r.listingOnly) { printListing(out, r); return; }
    for (auto& s : r.info) out << s << "\n";
    for (auto& lig : r.ligands) {
        for (auto& s : lig.info) out << s << "\n";
        out << fmt("\n=== Ligand %s %d chain %s: %zu atoms, %zu rings ===\n", lig.resName.c_str(), lig.resi,
                   lig.chain.c_str(), lig.atoms.size(), lig.nRings);
        out << fmt("%-24s %-3s %-20s %-20s %-6s %s\n", "Type", "Fav", "Ligand atom", "Protein atom", "d(A)", "Notes");
        std::map<std::string, int> cnt;
        for (auto& h : lig.hits) {
            out << fmt("%-24s %-3s %-20s %-20s %-6.2f %s\n", h.type.c_str(), h.favorable ? "+" : "-", h.ligandAtom.c_str(),
                       h.proteinAtom.c_str(), h.distance, h.note.c_str());
            cnt[h.type]++;
        }
        out << "\nSummary by type:\n";
        for (auto& kv : cnt) out << fmt("  %-28s %d\n", kv.first.c_str(), kv.second);

        std::map<std::string, std::pair<int, std::set<std::string>>> byRes;
        for (auto& h : lig.hits) { byRes[h.residue].first++; byRes[h.residue].second.insert(h.type); }
        std::vector<std::pair<std::string, std::pair<int, std::set<std::string>>>> rv(byRes.begin(), byRes.end());
        std::sort(rv.begin(), rv.end(), [](auto& a, auto& b) { return a.second.first > b.second.first; });
        out << "\nSummary by protein residue:\n";
        for (auto& e : rv) {
            std::string ts;
            for (auto& t : e.second.second) ts += (ts.empty() ? "" : ", ") + t;
            out << fmt("  %-12s %2d  %s\n", e.first.c_str(), e.second.first, ts.c_str());
        }
        if (!lig.metalSphere.empty())
            out << "\nCoordination sphere of the metals in contact with the ligand:\n" << lig.metalSphere;
    }
}

void printCsv(std::ostream& out, std::ostream& info, const Result& r) {
    if (r.listingOnly) { printListing(out, r); return; }
    out << "ligand,type,favorable,ligand_atom,protein_atom,distance,notes\n";
    for (auto& s : r.info) info << s << "\n";
    for (auto& lig : r.ligands) {
        for (auto& s : lig.info) info << s << "\n";
        for (auto& h : lig.hits)
            out << fmt("%s,%s,%d,%s,%s,%.2f,%s\n", lig.id.c_str(), h.type.c_str(), h.favorable ? 1 : 0, h.ligandAtom.c_str(),
                       h.proteinAtom.c_str(), h.distance, h.note.c_str());
    }
}

}  // namespace ligint
