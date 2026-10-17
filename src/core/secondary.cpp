// Protein backbone for the graphics and secondary-structure assignment (Kabsch & Sander, Biopolymers 22 (1983) 2577).
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "geom.hpp"
#include "internal.hpp"
#include "ligint/scene3d.hpp"

namespace ligint {
namespace {
using namespace geom;

// Hydrogen bond energy (kcal/mol) between the C=O of the acceptor residue a and the N-H of the donor residue b.
double hbondEnergy(const V3& c, const V3& o, const V3& n, const V3& h) {
    const double rON = dist(o, n), rCH = dist(c, h), rOH = dist(o, h), rCN = dist(c, n);
    if (rON < 0.5 || rCH < 0.5 || rOH < 0.5 || rCN < 0.5) return 0;
    return 27.888 * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN);
}

}  // namespace

void assignSecondaryStructure(std::vector<BackboneResidue>& r) {
    const int N = static_cast<int>(r.size());
    for (auto& x : r) x.ss = 'C';
    if (N < 4) return;

    // positions of the amide hydrogens: 1 A from N, opposite to the C=O of the previous residue
    std::vector<char> hasH(N, 0), complete(N, 0);
    std::vector<V3> H(N, V3{0, 0, 0});
    for (int i = 0; i < N; i++) {
        complete[i] = r[i].hasN && r[i].hasC && r[i].hasO;
        if (i > 0 && r[i].linkedToPrevious && complete[i] && r[i - 1].hasC && r[i - 1].hasO && r[i].res != "PRO") {
            const V3 d = normalized(sub(r[i - 1].c, r[i - 1].o));
            H[i] = add(r[i].n, d);
            hasH[i] = 1;
        }
    }

    // hydrogen bonds: key a*N+b  <=>  C=O of a accepts from N-H of b (energy < -0.5 kcal/mol)
    std::unordered_set<std::uint64_t> bonds;
    auto key = [&](int a, int b) { return static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(N) + static_cast<std::uint64_t>(b); };
    {
        const double cell = 9.0;
        std::unordered_map<std::int64_t, std::vector<int>> grid;
        auto cellKey = [&](int i, int j, int k) { return ((std::int64_t)(i + 4096) << 42) | ((std::int64_t)(j + 4096) << 21) | (std::int64_t)(k + 4096); };
        for (int i = 0; i < N; i++)
            grid[cellKey((int)std::floor(r[i].ca[0] / cell), (int)std::floor(r[i].ca[1] / cell), (int)std::floor(r[i].ca[2] / cell))].push_back(i);
        for (int a = 0; a < N; a++) {
            if (!complete[a]) continue;
            const int ci = (int)std::floor(r[a].ca[0] / cell), cj = (int)std::floor(r[a].ca[1] / cell), ck = (int)std::floor(r[a].ca[2] / cell);
            for (int i = ci - 1; i <= ci + 1; i++)
                for (int j = cj - 1; j <= cj + 1; j++)
                    for (int k = ck - 1; k <= ck + 1; k++) {
                        auto it = grid.find(cellKey(i, j, k));
                        if (it == grid.end()) continue;
                        for (int b : it->second) {
                            if (b == a || !hasH[b] || std::abs(a - b) < 2) continue;
                            if (dist(r[a].ca, r[b].ca) > 9.0) continue;
                            if (hbondEnergy(r[a].c, r[a].o, r[b].n, H[b]) < -0.5) bonds.insert(key(a, b));
                        }
                    }
        }
    }
    auto hb = [&](int a, int b) { return a >= 0 && b >= 0 && a < N && b < N && bonds.count(key(a, b)) > 0; };
    // the residues i..i+n are one uninterrupted stretch of chain
    auto continuous = [&](int i, int n) {
        if (i < 0 || i + n >= N) return false;
        for (int k = 1; k <= n; k++) if (!r[i + k].linkedToPrevious) return false;
        return true;
    };

    // helices: two consecutive n-turns at i-1 and i make the residues i..i+n-1 a helix (n = 3: 3-10, 4: alpha, 5: pi)
    std::vector<char> helix(N, 0);
    for (int n = 3; n <= 5; n++) {
        std::vector<char> turn(N, 0);
        for (int i = 0; i + n < N; i++) turn[i] = continuous(i, n) && hb(i, i + n);
        for (int i = 1; i + n - 1 < N; i++)
            if (turn[i - 1] && turn[i]) for (int k = 0; k < n; k++) helix[i + k] = 1;
    }

    // bridges between non-adjacent residues; consecutive bridged residues form strands
    std::vector<char> bridge(N, 0);
    {
        const double cell = 9.0;
        std::unordered_map<std::int64_t, std::vector<int>> grid;
        auto cellKey = [&](int i, int j, int k) { return ((std::int64_t)(i + 4096) << 42) | ((std::int64_t)(j + 4096) << 21) | (std::int64_t)(k + 4096); };
        for (int i = 0; i < N; i++)
            grid[cellKey((int)std::floor(r[i].ca[0] / cell), (int)std::floor(r[i].ca[1] / cell), (int)std::floor(r[i].ca[2] / cell))].push_back(i);
        for (int i = 1; i + 1 < N; i++) {
            if (!continuous(i - 1, 2)) continue;
            const int ci = (int)std::floor(r[i].ca[0] / cell), cj = (int)std::floor(r[i].ca[1] / cell), ck = (int)std::floor(r[i].ca[2] / cell);
            for (int a = ci - 1; a <= ci + 1; a++)
                for (int b = cj - 1; b <= cj + 1; b++)
                    for (int c = ck - 1; c <= ck + 1; c++) {
                        auto it = grid.find(cellKey(a, b, c));
                        if (it == grid.end()) continue;
                        for (int j : it->second) {
                            if (j <= i + 2 || j + 1 >= N || !continuous(j - 1, 2)) continue;
                            const bool parallel = (hb(i - 1, j) && hb(j, i + 1)) || (hb(j - 1, i) && hb(i, j + 1));
                            const bool anti = (hb(i, j) && hb(j, i)) || (hb(i - 1, j + 1) && hb(j - 1, i + 1));
                            if (parallel || anti) bridge[i] = bridge[j] = 1;
                        }
                    }
        }
    }
    std::vector<char> strand(N, 0);
    for (int i = 0; i < N; i++) {
        if (!bridge[i]) continue;
        const bool prev = i > 0 && r[i].linkedToPrevious && bridge[i - 1];
        const bool next = i + 1 < N && r[i + 1].linkedToPrevious && bridge[i + 1];
        if (prev || next) strand[i] = 1;
    }
    for (int i = 1; i + 1 < N; i++)     // a single residue between two strand residues (beta bulge) belongs to the strand
        if (!strand[i] && strand[i - 1] && strand[i + 1] && r[i].linkedToPrevious && r[i + 1].linkedToPrevious) strand[i] = 1;

    for (int i = 0; i < N; i++) r[i].ss = helix[i] ? 'H' : (strand[i] ? 'E' : 'C');
}

namespace detail {

// Builds Result::backbone (with the secondary structure) and Result::proteinAtoms from the protein atoms (file order).
void buildProteinScene(const std::vector<Atom>& prot, Result& res) {
    struct Raw {
        std::string res, chain;
        int resi = 0;
        char icode = ' ';
        V3 n{0, 0, 0}, ca{0, 0, 0}, c{0, 0, 0}, o{0, 0, 0};
        bool hn = false, hca = false, hc = false, ho = false, het = true;
        double bf = 0;
    };
    std::vector<Raw> list;
    std::map<std::tuple<std::string, int, char>, int> where;
    std::vector<int> atomRaw(prot.size(), -1);
    for (size_t i = 0; i < prot.size(); i++) {
        const Atom& a = prot[i];
        const auto key = std::make_tuple(a.chain, a.resi, a.icode);
        auto it = where.find(key);
        int idx;
        if (it == where.end()) {
            idx = static_cast<int>(list.size());
            where[key] = idx;
            Raw r; r.res = a.res; r.chain = a.chain; r.resi = a.resi; r.icode = a.icode;
            list.push_back(r);
        } else {
            idx = it->second;
        }
        atomRaw[i] = idx;
        Raw& r = list[idx];
        if (!a.het) r.het = false;
        if (a.name == "N") { r.n = detail::pos(a); r.hn = true; }
        else if (a.name == "CA" && a.el == "C") { r.ca = detail::pos(a); r.hca = true; r.bf = a.bf; }
        else if (a.name == "C") { r.c = detail::pos(a); r.hc = true; }
        else if (a.name == "O") { r.o = detail::pos(a); r.ho = true; }
    }

    // a residue belongs to the backbone with its C-alpha; a HETATM residue (selenomethionine...) needs the whole peptide unit
    std::vector<int> toBackbone(list.size(), -1);
    std::map<std::string, int> chains;
    res.backbone.clear();
    for (size_t i = 0; i < list.size(); i++) {
        const Raw& r = list[i];
        if (!r.hca || (r.het && !(r.hn && r.hc && r.ho))) continue;
        BackboneResidue b;
        b.res = r.res; b.chain = r.chain; b.resi = r.resi; b.icode = r.icode;
        b.n = r.n; b.ca = r.ca; b.c = r.c; b.o = r.o;
        b.hasN = r.hn; b.hasC = r.hc; b.hasO = r.ho;
        b.bfactor = r.bf;
        auto ch = chains.find(r.chain);
        if (ch == chains.end()) ch = chains.emplace(r.chain, static_cast<int>(chains.size())).first;
        b.chainIndex = ch->second;
        if (!res.backbone.empty()) {
            const BackboneResidue& p = res.backbone.back();
            if (p.chain == b.chain) {
                if (p.hasC && b.hasN) b.linkedToPrevious = geom::dist(p.c, b.n) < 2.0;
                else b.linkedToPrevious = geom::dist(p.ca, b.ca) < 4.3;
            }
        }
        toBackbone[i] = static_cast<int>(res.backbone.size());
        res.backbone.push_back(b);
    }
    assignSecondaryStructure(res.backbone);

    res.proteinAtoms.clear();
    res.proteinAtoms.reserve(prot.size());
    for (size_t i = 0; i < prot.size(); i++) {
        const Atom& a = prot[i];
        const int bb = toBackbone[atomRaw[i]];
        ProteinAtom p;
        p.pos = detail::pos(a);
        p.radius = static_cast<float>(detail::vdw(a.el));
        p.element[0] = a.el.size() > 0 ? a.el[0] : 'C';
        p.element[1] = a.el.size() > 1 ? a.el[1] : 0;
        p.element[2] = 0;
        p.residue = bb;
        p.hetero = a.het && bb < 0;
        res.proteinAtoms.push_back(p);
    }
}

}  // namespace detail
}  // namespace ligint
