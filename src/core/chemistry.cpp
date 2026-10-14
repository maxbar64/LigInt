// Rings and chemical roles (donor, acceptor, charge...) of ligand and protein.
#include "internal.hpp"

namespace ligint {
namespace detail {

// ---------- rings ----------
bool makeRing(const vector<const Atom*>& pts, const string& lab, const string& pres,
                     double tol, Ring& out) {
    size_t m = pts.size();
    V3 c = {0, 0, 0}, n = {0, 0, 0};
    for (auto p : pts) { V3 q = pos(*p); for (int i = 0; i < 3; i++) c[i] += q[i] / m; }
    for (size_t i = 0; i < m; i++) {  // Newell's method
        V3 a = pos(*pts[i]), b = pos(*pts[(i + 1) % m]);
        n[0] += (a[1]-b[1]) * (a[2]+b[2]);
        n[1] += (a[2]-b[2]) * (a[0]+b[0]);
        n[2] += (a[0]-b[0]) * (a[1]+b[1]);
    }
    double l = nrm(n);
    if (l < 1e-6) return false;
    for (int i = 0; i < 3; i++) n[i] /= l;
    for (auto p : pts) if (fabs(dot(sub(pos(*p), c), n)) > tol) return false;
    out = {c, n, lab, pres};
    return true;
}
vector<vector<int>> findRings(const vector<Atom>& L) {  // cycles of 5-6 atoms
    set<vector<int>> seen; vector<vector<int>> out;
    function<void(int,int,vector<int>&)> dfs = [&](int s, int cur, vector<int>& path) {
        for (int nx : L[cur].nb) {
            if (isMetal(L[nx].el)) continue;
            if (nx == s && path.size() >= 5) {
                auto key = path; sort(key.begin(), key.end());
                if (seen.insert(key).second) out.push_back(path);
            } else if (nx > s && path.size() < 6 && find(path.begin(), path.end(), nx) == path.end()) {
                path.push_back(nx); dfs(s, nx, path); path.pop_back();
            }
        }
    };
    for (size_t s = 0; s < L.size(); s++) { vector<int> p = {(int)s}; dfs(s, s, p); }
    return out;
}

// ---------- roles ----------
bool isNuc(const string& R) {
    static const set<string> s = {"A","C","G","U","T","I","DA","DC","DG","DT","DU","DI"};
    return s.count(R);
}
Role proteinRole(const Atom& a) {
    Role r;
    const string &R = a.res, &N = a.name;
    if (isMetal(a.el)) { r.metal = true; return r; }
    if (isNuc(R)) {  // DNA / RNA: phosphate, sugar, bases
        string n = N;
        for (auto& c : n) if (c == '*') c = '\'';
        char b = R.back();
        if (n == "OP1" || n == "OP2" || n == "O1P" || n == "O2P") { r.acc = true; r.q = -1; }
        else if (n == "O5'" || n == "O3'" || n == "O4'") r.acc = true;
        else if (n == "O2'") r.don = r.acc = true;
        else if (b == 'A') { if (n == "N6") r.don = true; if (n == "N1" || n == "N3" || n == "N7") r.acc = true; }
        else if (b == 'G' || b == 'I') {
            if (n == "N1" || (b == 'G' && n == "N2")) r.don = true;
            if (n == "O6" || n == "N3" || n == "N7") r.acc = true;
        }
        else if (b == 'C') { if (n == "N4") r.don = true; if (n == "O2" || n == "N3") r.acc = true; }
        else if (b == 'T' || b == 'U') { if (n == "N3") r.don = true; if (n == "O2" || n == "O4") r.acc = true; }
        if (a.el == "C" && (n == "C7" || n == "C5M" || (n == "C5" && (b == 'T' || b == 'U' || b == 'C')))) r.hyd = true;
        return r;
    }
    if (N == "N" && R != "PRO") r.don = true;
    if (N == "O" || N == "OXT") r.acc = true;
    auto in = [&](initializer_list<const char*> l) { for (auto s : l) if (N == s) return true; return false; };
    if ((R=="SER" && N=="OG") || (R=="THR" && N=="OG1") || (R=="TYR" && N=="OH")) r.don = r.acc = true;
    if (R=="ASN") { if (N=="OD1") r.acc = true; if (N=="ND2") r.don = true; }
    if (R=="GLN") { if (N=="OE1") r.acc = true; if (N=="NE2") r.don = true; }
    if (R=="ASP" && in({"OD1","OD2"})) { r.acc = true; r.q = -1; }
    if (R=="GLU" && in({"OE1","OE2"})) { r.acc = true; r.q = -1; }
    if (R=="LYS" && N=="NZ") { r.don = true; r.q = +1; }
    if (R=="ARG" && in({"NE","NH1","NH2"})) { r.don = true; r.q = +1; }
    if (R=="HIS" && in({"ND1","NE2"})) r.don = r.acc = true;
    if (R=="TRP" && N=="NE1") r.don = true;
    if (a.el == "C" && N != "C" && N != "CA") {
        if (R=="ALA"||R=="VAL"||R=="LEU"||R=="ILE"||R=="MET"||R=="PRO"||R=="PHE"||R=="TRP"||R=="CYS") r.hyd = true;
        if (R=="TYR" && in({"CB","CG","CD1","CD2","CE1","CE2"})) r.hyd = true;
        if (R=="LYS" && in({"CB","CG","CD"})) r.hyd = true;
        if (R=="ARG" && in({"CB","CG"})) r.hyd = true;
        if (R=="HIS" && in({"CB","CG","CD2"})) r.hyd = true;
        if (R=="THR" && N=="CG2") r.hyd = true;
    }
    // with explicit hydrogens in the PDB (applied by the analysis) the roles are refined
    return r;
}
void refineProteinRoleWithH(const Atom& a, Role& r) {
    if (r.metal || !isNO(a.el)) return;
    bool hasH = !a.hs.empty();
    r.don = hasH;
    if (a.res == "HIS" && hasH && (a.name == "ND1" || a.name == "NE2")) r.acc = false;
    if (hasH && (a.res == "ASP" || a.res == "GLU") && r.q < 0) r.q = 0;   // protonated carboxyl
    if (!hasH && a.res == "LYS" && a.name == "NZ") r.q = 0;               // neutral lysine
    if (a.res == "ARG" && !hasH) r.q = 0;
}

// is there a path u -> v of at most maxD bonds that does not go through 'avoid'? (detects ring membership)
static bool connectedAvoiding(const vector<Atom>& L, int u, int v, int avoid, int maxD) {
    vector<int> dep(L.size(), -1);
    deque<int> q;
    dep[u] = 0; q.push_back(u);
    while (!q.empty()) {
        int c = q.front(); q.pop_front();
        if (c == v) return true;
        if (dep[c] >= maxD) continue;
        for (int n : L[c].nb) if (n != avoid && dep[n] < 0) { dep[n] = dep[c] + 1; q.push_back(n); }
    }
    return false;
}
// sp3 carbon (all bonds as long as single bonds); used to recognize aliphatic amines
static bool sp3C(const vector<Atom>& L, int k) {
    const Atom& c = L[k];
    if (c.el != "C") return false;
    for (int j : c.nb) {
        double d = dist(c, L[j]);
        const string& e = L[j].el;
        double lim = e == "C" ? 1.48 : e == "N" ? 1.41 : e == "O" ? 1.40 : 1.60;
        if (d < lim) return false;
    }
    return true;
}
static bool inRingC(const vector<Atom>& L, int k) {  // is atom k part of a ring (cycles up to 7 atoms)?
    const Atom& c = L[k];
    for (size_t x = 0; x < c.nb.size(); x++)
        for (size_t y = x + 1; y < c.nb.size(); y++)
            if (connectedAvoiding(L, c.nb[x], c.nb[y], k, 5)) return true;
    return false;
}
static bool termO(const vector<Atom>& L, int k) { return L[k].el == "O" && L[k].nb.size() == 1; }

Role ligandRole(const vector<Atom>& L, int i, bool hKnown) {
    const Atom& a = L[i];
    Role r;
    size_t deg = a.nb.size(), nh = a.hs.size();
    if (isMetal(a.el)) { r.metal = true; return r; }
    if (a.el == "Cl" || a.el == "Br" || a.el == "I") r.hal = true;
    if (a.fc != 0) r.q = a.fc > 0 ? +1 : -1;

    if (a.el == "O") {
        if (deg == 0) r.don = r.acc = true;
        else if (deg >= 2) r.acc = true;          // ether / ester / phosphate bridge: acceptor only
        else {
            const Atom& X = L[a.nb[0]];
            r.acc = true;
            if (hKnown) r.don = nh > 0;
            else {
                double d = dist(a, X);
                r.don = (X.el == "C" || X.el == "N") ? d > 1.32 : d > 1.52;   // C-OH / N-OH / P-OH versus C=O / N=O / P=O
            }
            if (X.el == "C" || X.el == "P" || X.el == "S") {
                vector<int> sib;
                for (int k : X.nb) if (termO(L, k)) sib.push_back(k);
                bool anion = sib.size() >= 2;      // carboxylate / phosphate / sulfonate (physiological pH)
                if (anion && hKnown) {
                    if (X.el == "C") { for (int k : sib) if (!L[k].hs.empty()) anion = false; }
                    else if (nh > 0) anion = false;
                }
                if (anion) { r.q = -1; r.don = false; }
            }
        }
    }
    if (a.el == "N") {
        int tO = 0;
        for (int k : a.nb) if (termO(L, k)) tO++;
        if (tO >= 2) {                        // nitro: neutral, weak acceptor -> ignored
            r.q = 0;
        } else if (hKnown) {
            r.don = nh > 0;
            size_t tot = deg + nh;
            if (a.fc > 0 || tot >= 4) { r.q = +1; r.acc = false; }
            else if (tot <= 2) r.acc = true;                 // sp2 N with a lone pair / nitrile
            else {                                           // 3 substituents: planar (amide, pyrrole) or pyramidal
                vector<V3> s;
                for (int k : a.nb) s.push_back(pos(L[k]));
                for (auto& h : a.hs) s.push_back(h);
                double sum = angV(s[0], pos(a), s[1]) + angV(s[0], pos(a), s[2]) + angV(s[1], pos(a), s[2]);
                r.acc = sum < 350.0;
            }
        } else {
            if (deg <= 2) r.don = r.acc = true;
            else if (deg == 3) {
                double s = angleDeg(L[a.nb[0]], a, L[a.nb[1]]) + angleDeg(L[a.nb[0]], a, L[a.nb[2]]) +
                           angleDeg(L[a.nb[1]], a, L[a.nb[2]]);
                if (s < 345.0) r.acc = true;
            } else { r.don = true; r.q = +1; }
        }
        if (!hKnown) {  // N bonded to a carbonyl C (amide, carbamate, urea): the lone pair is delocalized
            for (int k : a.nb) {
                const Atom& C = L[k];
                if (C.el != "C") continue;
                for (int j : C.nb)
                    if (L[j].el == "O" && L[j].nb.size() == 1 && dist(C, L[j]) < 1.30) r.acc = false;
            }
        }
        if (!hKnown && tO < 2 && deg >= 1 && deg <= 3) {  // aliphatic amine: protonated at physiological pH
            bool amine = true;
            for (int k : a.nb) if (!sp3C(L, k)) { amine = false; break; }
            if (amine) { r.q = +1; r.don = true; r.acc = false; }
        }
        for (int k : a.nb) {                  // guanidinium / amidinium / imidazolium: delocalized charge
            const Atom& C = L[k];
            if (C.el != "C" || C.nb.size() < 2) continue;
            int nN = 0; bool anyCat = false;
            for (int j : C.nb) if (L[j].el == "N") { nN++; if (L[j].fc > 0) anyCat = true; }
            if (nN >= 2 && anyCat && (nh > 0 || a.fc > 0)) { r.q = +1; r.acc = false; if (nh > 0) r.don = true; }
            if (!hKnown && C.nb.size() == 3 && nN == 3 && !inRingC(L, k)) { r.q = +1; r.don = true; r.acc = false; }  // guanidinium (not in a ring)
            if (!hKnown && C.nb.size() == 3 && nN == 2) {  // amidine R-C(N)=N: protonated (amidinium) at physiological pH
                int third = -1; vector<int> ns;
                for (int j : C.nb) { if (L[j].el == "N") ns.push_back(j); else third = j; }
                bool amid = third >= 0 && L[third].el == "C";
                for (int j : ns) if (dist(C, L[j]) > 1.38 || L[j].nb.size() > 2) amid = false;  // short C-N bonds, N not substituted twice
                if (amid) {                                                                      // C not in a ring (excludes aminopyridines, imidazoles...)
                    int ci = k;
                    for (size_t x = 0; x < C.nb.size() && amid; x++)
                        for (size_t y = x + 1; y < C.nb.size() && amid; y++)
                            if (connectedAvoiding(L, C.nb[x], C.nb[y], ci, 5)) amid = false;
                }
                if (amid) { r.q = +1; r.don = true; r.acc = false; }
            }
        }
    }
    if (a.el == "C" && a.fc == 0) {
        bool h = true;
        for (int k : a.nb) if (isNO(L[k].el) || L[k].el == "P") h = false;
        r.hyd = h;
    }
    return r;
}

}  // namespace detail
}  // namespace ligint
