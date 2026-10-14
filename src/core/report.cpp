// Rapporti testuali: tabella leggibile e CSV (stesso formato del programma originale).
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
    out << "Ligandi (HETATM) presenti:";
    for (auto& kv : r.hetCounts) out << " " << kv.first << "(" << kv.second << ")";
    out << "\n";
    if (!r.notFound.empty()) out << "Nessun ligando '" << r.notFound << "' trovato.\n";
}

}  // namespace

void printText(std::ostream& out, const Result& r) {
    if (r.listingOnly) { printListing(out, r); return; }
    for (auto& s : r.info) out << s << "\n";
    for (auto& lig : r.ligands) {
        for (auto& s : lig.info) out << s << "\n";
        out << fmt("\n=== Ligando %s %d catena %s: %zu atomi, %zu anelli ===\n", lig.resName.c_str(), lig.resi,
                   lig.chain.c_str(), lig.atoms.size(), lig.nRings);
        out << fmt("%-24s %-3s %-20s %-20s %-6s %s\n", "Tipo", "F/S", "Atomo ligando", "Atomo proteina", "d(A)", "Note");
        std::map<std::string, int> cnt;
        for (auto& h : lig.hits) {
            out << fmt("%-24s %-3s %-20s %-20s %-6.2f %s\n", h.type.c_str(), h.favorable ? "+" : "-", h.ligandAtom.c_str(),
                       h.proteinAtom.c_str(), h.distance, h.note.c_str());
            cnt[h.type]++;
        }
        out << "\nRiepilogo per tipo:\n";
        for (auto& kv : cnt) out << fmt("  %-28s %d\n", kv.first.c_str(), kv.second);

        std::map<std::string, std::pair<int, std::set<std::string>>> byRes;
        for (auto& h : lig.hits) { byRes[h.residue].first++; byRes[h.residue].second.insert(h.type); }
        std::vector<std::pair<std::string, std::pair<int, std::set<std::string>>>> rv(byRes.begin(), byRes.end());
        std::sort(rv.begin(), rv.end(), [](auto& a, auto& b) { return a.second.first > b.second.first; });
        out << "\nRiepilogo per residuo di proteina:\n";
        for (auto& e : rv) {
            std::string ts;
            for (auto& t : e.second.second) ts += (ts.empty() ? "" : ", ") + t;
            out << fmt("  %-12s %2d  %s\n", e.first.c_str(), e.second.first, ts.c_str());
        }
        if (!lig.metalSphere.empty())
            out << "\nSfera di coordinazione dei metalli a contatto con il ligando:\n" << lig.metalSphere;
    }
}

void printCsv(std::ostream& out, std::ostream& info, const Result& r) {
    if (r.listingOnly) { printListing(out, r); return; }
    out << "ligando,tipo,favorevole,atomo_ligando,atomo_proteina,distanza,note\n";
    for (auto& s : r.info) info << s << "\n";
    for (auto& lig : r.ligands) {
        for (auto& s : lig.info) info << s << "\n";
        for (auto& h : lig.hits)
            out << fmt("%s,%s,%d,%s,%s,%.2f,%s\n", lig.id.c_str(), h.type.c_str(), h.favorable ? 1 : 0, h.ligandAtom.c_str(),
                       h.proteinAtom.c_str(), h.distance, h.note.c_str());
    }
}

}  // namespace ligint
