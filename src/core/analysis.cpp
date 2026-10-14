// Analisi: lettura del PDB, ruoli chimici, ricerca delle interazioni.
#include "internal.hpp"

namespace ligint {

// ---------------------------------------------------------------------------------------------
// Params
// ---------------------------------------------------------------------------------------------
bool Params::set(const std::string& k, double v) {
    if (k == "hb-min") hbMin = v; else if (k == "hb-max") hbMax = v;
    else if (k == "hb-ha") hbHA = v; else if (k == "hb-angle") hbDHA = v;
    else if (k == "hyd-min") hydMin = v; else if (k == "hyd-max") hydMax = v;
    else if (k == "salt-max") saltMax = v; else if (k == "rep-max") repMax = v;
    else if (k == "hal-max") halMax = v; else if (k == "hal-angle") halAngle = v;
    else if (k == "metal-max") metalMax = v; else if (k == "pi-max") piMax = v;
    else if (k == "pi-offset") piOffset = v; else if (k == "pi-angle") piAngle = v; else if (k == "cp-max") cpMax = v;
    else if (k == "clash") clash = v;
    else if (k == "hb-long") hbLong = v; else if (k == "hb-long-tol") hbLongTol = v;
    else if (k == "salt-centroid") saltCentroid = v; else if (k == "hb-ante") hbAnte = v;
    else if (k == "wb-min") wbMin = v; else if (k == "wb-max") wbMax = v;
    else return false;
    return true;
}

void Params::applyPlipLike() {
    hbLong = 4.1;
    saltCentroid = 5.5;
    halMax = 4.0;
}

std::vector<std::string> Params::clampToGrid() {
    std::vector<std::string> warnings;
    auto clampv = [&](double& x, const char* name) {
        if (x > kMaxGridRadius) {
            char b[160];
            std::snprintf(b, sizeof b, "Attenzione: --%s %g limitato a %g A (massimo della ricerca a griglia)", name, x, kMaxGridRadius);
            warnings.push_back(b);
            x = kMaxGridRadius;
        }
    };
    clampv(hbMax, "hb-max"); clampv(hbLong, "hb-long"); clampv(hydMax, "hyd-max"); clampv(saltMax, "salt-max");
    clampv(repMax, "rep-max"); clampv(halMax, "hal-max"); clampv(metalMax, "metal-max");
    return warnings;
}

// ---------------------------------------------------------------------------------------------
// Analisi
// ---------------------------------------------------------------------------------------------
namespace {
using namespace detail;
using namespace std;
typedef tuple<string,int,char> LKey;

AtomRecord toRecord(const Atom& a) {
    AtomRecord r;
    r.name = a.name; r.res = a.res; r.chain = a.chain; r.element = a.el;
    r.icode = a.icode; r.resi = a.resi; r.pos = pos(a); r.hetero = a.het;
    return r;
}
}  // namespace

Result analyzeFile(const std::string& path, const Options& opt) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Impossibile aprire " + path);
    return analyzeStream(in, opt);
}

Result analyzeStream(std::istream& in, const Options& opt) {
    using namespace detail;
    using namespace std;
    Result res;
    Params P = opt.params;
    for (auto& w : P.clampToGrid()) res.warnings.push_back(w);
    const bool useWater = opt.useWater, metalDonors = opt.metalDonors, backboneHyd = opt.backboneHyd, mergeLinked = opt.mergeLinked;
    const vector<Tmpl>& templates = opt.templates;
    const string& ligName = opt.ligName;
    const string& wantChain = opt.chain;
    const int wantResi = opt.resi;

    vector<Atom> prot, waters;
    vector<pair<Atom,bool>> hatoms;  // idrogeni del PDB (bool = appartiene a un ligando bersaglio)
    map<LKey, vector<Atom>> ligs;
    map<LKey, vector<Atom>> hetRes;   // tutti i residui HETATM (acqua esclusa), per unire quelli legati
    map<string, set<LKey>> hetNames;
    map<string, char> altSeen;
    string line;
    long nBadCoord = 0;
    while (getline(in, line)) {
        if (line.compare(0, 6, "ENDMDL") == 0) break;  // solo primo modello
        if (line.size() < 54) continue;
        string rec = line.substr(0, 6);
        if (rec != "ATOM  " && rec != "HETATM") continue;
        Atom a;
        a.name = trim(line.substr(12, 4));
        a.res = trim(line.substr(17, 3));
        a.chain = string(1, line[21]);
        a.resi = atoi(line.substr(22, 4).c_str());
        a.icode = line[26];
        char alt = line[16];
        if (alt != ' ') {  // per ogni residuo si tiene la prima conformazione incontrata
            string rk = line.substr(21, 6) + line.substr(17, 3);
            auto it = altSeen.find(rk);
            if (it == altSeen.end()) altSeen[rk] = alt;
            else if (it->second != alt) continue;
        }
        a.x = atof(line.substr(30, 8).c_str());
        a.y = atof(line.substr(38, 8).c_str());
        a.z = atof(line.substr(46, 8).c_str());
        if (!isfinite(a.x) || !isfinite(a.y) || !isfinite(a.z) || fabs(a.x) > 1e5 || fabs(a.y) > 1e5 || fabs(a.z) > 1e5) {
            nBadCoord++;
            continue;
        }
        a.bf = line.size() >= 66 ? atof(line.substr(60, 6).c_str()) : 0.0;
        string el = line.size() >= 78 ? trim(line.substr(76, 2)) : "";
        if (el.empty()) {
            string nm = line.substr(12, 4);
            if (rec == "HETATM" && a.name.size() <= 2 && isalpha((unsigned char)nm[0]) && isalpha((unsigned char)nm[1]))
                el = nm.substr(0, 2);
            else for (char c : a.name) if (isalpha((unsigned char)c)) { el = string(1, c); break; }
        }
        if (el.empty()) continue;
        el = normEl(el);
        a.el = el;
        a.het = rec == "HETATM";
        bool isH = isHel(el);
        if (a.res == "HOH" || a.res == "WAT" || a.res == "DOD") {
            if (useWater && el == "O") waters.push_back(a);
            continue;
        }
        bool isLig = rec == "HETATM" && a.res == ligName &&
                     (wantChain.empty() || a.chain == wantChain) &&
                     (wantResi == INT_MIN || a.resi == wantResi);
        if (isH) { hatoms.push_back({a, isLig}); continue; }
        if (rec == "HETATM") { hetNames[a.res].insert({a.chain, a.resi, a.icode}); hetRes[{a.chain, a.resi, a.icode}].push_back(a); }
        if (isLig) ligs[{a.chain, a.resi, a.icode}].push_back(a);
        else prot.push_back(a);
    }
    if (nBadCoord > 0) res.warnings.push_back("Attenzione: " + to_string(nBadCoord) + " atomi con coordinate non valide ignorati");
    if (ligName.empty() || ligs.empty()) {
        res.listingOnly = true;
        for (auto& [n, s] : hetNames) res.hetCounts[n] = (int)s.size();
        for (auto& [k, atoms] : hetRes) {
            LigandId id;
            id.res = atoms[0].res; id.chain = get<0>(k); id.resi = get<1>(k); id.icode = get<2>(k); id.atoms = (int)atoms.size();
            res.hetResidues.push_back(id);
        }
        if (!ligName.empty()) res.notFound = ligName;
        return res;
    }

    Grid grid; grid.build(prot);
    vector<int> tmp;

    // idrogeni del PDB: proteina -> atomo pesante piu' vicino nello stesso residuo; ligando -> da parte
    map<LKey, vector<V3>> ligH;
    int nHprot = 0;
    for (auto& [h, isl] : hatoms) {
        if (isl) { ligH[{h.chain, h.resi, h.icode}].push_back(pos(h)); continue; }
        grid.near(pos(h), 1.3, tmp);
        int best = -1; double bd = 1e9;
        for (int k : tmp) {
            const Atom& q = prot[k];
            if (q.chain != h.chain || q.resi != h.resi || q.icode != h.icode) continue;
            double d = nrm(sub(pos(q), pos(h)));
            if (d < bd) { bd = d; best = k; }
        }
        if (best >= 0) { prot[best].hs.push_back(pos(h)); nHprot++; }
    }
    bool hKnownProt = nHprot >= 10;

    vector<Role> pr(prot.size());
    for (size_t i = 0; i < prot.size(); i++) {
        pr[i] = proteinRole(prot[i]);
        if (hKnownProt) refineProteinRoleWithH(prot[i], pr[i]);
    }

    // N/O della proteina che coordinano un metallo (es. His di uno zinc finger) non hanno H da donare
    if (!hKnownProt) {
        for (size_t k = 0; k < prot.size(); k++) {
            if (!isNO(prot[k].el) || pr[k].metal || (!pr[k].don && !pr[k].acc)) continue;
            grid.near(pos(prot[k]), 3.0, tmp);
            for (int m : tmp)
                if (pr[m].metal && dist(prot[m], prot[k]) <= metalLimit(prot[m].el, prot[k].el)) {
                    pr[k].don = false;
                    if (prot[k].el == "N") pr[k].acc = false;   // His/N legata al metallo: coppia solitaria impegnata
                    break;
                }
        }
    }

    // connettivita' della proteina, calcolata su richiesta
    unordered_map<int, vector<int>> pnbCache;
    auto pnb = [&](int k) -> const vector<int>& {
        auto it = pnbCache.find(k);
        if (it != pnbCache.end()) return it->second;
        vector<int> o, c2;
        grid.near(pos(prot[k]), 2.2, c2);
        for (int j : c2) if (j != k && dist(prot[k], prot[j]) < covr(prot[k].el) + covr(prot[j].el) + 0.4) o.push_back(j);
        return pnbCache[k] = o;
    };
    auto angleOk = [&](const Atom& c, const vector<int>& nbs, const vector<Atom>& arr, const Atom& partner) {
        for (int j : nbs) if (angleDeg(arr[j], c, partner) < P.hbAnte) return false;
        return true;
    };
    struct HBres { bool ok = false; string note; };
    auto donorTest = [&](const Atom& D, const vector<int>& dnb, const vector<Atom>& arr, const Atom& A) {
        HBres r;
        if (!D.hs.empty()) {
            double bestAng = -1, bestD = 0;
            for (auto& h : D.hs) {
                double dha = nrm(sub(h, pos(A)));
                if (dha > P.hbHA) continue;
                double ang = angV(pos(D), h, pos(A));
                if (ang >= P.hbDHA && ang > bestAng) { bestAng = ang; bestD = dha; }
            }
            if (bestAng >= 0) {
                char b[80]; snprintf(b, sizeof b, "H...A %.2f; D-H...A %.0f gradi", bestD, bestAng);
                r.ok = true; r.note = b;
            }
        } else if (dnb.size() == 3) {  // ammina a tre sostituenti: l'H (coppia solitaria) punta opposto ai sostituenti
            V3 u = {0, 0, 0};
            for (int j : dnb) {
                V3 v = sub(pos(arr[j]), pos(D)); double n = nrm(v);
                for (int i = 0; i < 3; i++) u[i] += v[i] / n;
            }
            double nu = nrm(u);
            if (nu < 0.3) r.ok = angleOk(D, dnb, arr, A);   // planare (ammide, anilina): criterio generico
            else {
                V3 w = sub(pos(A), pos(D));
                double ang = acos(max(-1.0, min(1.0, -dot(u, w) / (nu * nrm(w))))) * 180.0 / M_PI;
                if (ang <= 50.0) { char b[80]; snprintf(b, sizeof b, "asse N-H atteso: %.0f gradi", ang); r.ok = true; r.note = b; }
            }
        } else r.ok = angleOk(D, dnb, arr, A);
        return r;
    };
    // H-bond LUNGO (D...A oltre hb-max): serve che l'H possa davvero puntare verso l'accettore
    auto donorTestLong = [&](const Atom& D, const vector<int>& dnb, const vector<Atom>& arr, const Atom& A) {
        HBres r;
        V3 w = sub(pos(A), pos(D));
        double nw = nrm(w);
        if (nw < 1e-6) return r;
        char b[80];
        if (!D.hs.empty()) {                                   // H noti: angolo D-H...A reale
            double best = -1;
            for (auto& h : D.hs) best = max(best, angV(pos(D), h, pos(A)));
            if (best >= 180.0 - P.hbLongTol) { snprintf(b, sizeof b, "lungo; D-H...A %.0f gradi", best); r.ok = true; r.note = b; }
        } else if (dnb.size() >= 2) {                          // H in posizione fissa, opposta ai sostituenti
            V3 u = {0, 0, 0};
            for (int j : dnb) {
                V3 v = sub(pos(arr[j]), pos(D)); double n = nrm(v);
                for (int i = 0; i < 3; i++) u[i] += v[i] / n;
            }
            double nu = nrm(u);
            if (nu > 0.3) {                                    // planare ammide N-H, ammina sp3...
                double ang = acos(max(-1.0, min(1.0, -dot(u, w) / (nu * nw)))) * 180.0 / M_PI;
                if (ang <= P.hbLongTol) { snprintf(b, sizeof b, "lungo; asse N-H atteso: %.0f gradi", ang); r.ok = true; r.note = b; }
            }
        } else if (dnb.size() == 1) {                          // idrossile, NH2, NH3+: l'H ruota attorno al legame
            const Atom& X = arr[dnb[0]];
            V3 v1 = sub(pos(X), pos(D));
            double th = acos(max(-1.0, min(1.0, dot(v1, w) / (nrm(v1) * nw)))) * 180.0 / M_PI;   // X-D...A
            double th0 = (D.el == "O" || dist(D, X) >= 1.42) ? 109.5 : 120.0;                      // angolo X-D-H atteso
            if (fabs(th - th0) <= P.hbLongTol) { snprintf(b, sizeof b, "lungo; X-D...A %.0f gradi (atteso %.0f)", th, th0); r.ok = true; r.note = b; }
        } else { r.ok = true; r.note = "lungo; donatore senza vicini"; }
        return r;
    };

    // gruppi carichi della proteina (centroide per residuo e segno), per --salt-centroid
    struct QGroup { V3 c{0, 0, 0}; int q = 0; vector<int> idx; };
    vector<QGroup> pQ;
    if (P.saltCentroid > 0) {
        map<tuple<string,int,char,string,int>, vector<int>> mq;
        for (size_t k = 0; k < prot.size(); k++)
            if (!prot[k].het && !pr[k].metal && pr[k].q != 0)
                mq[{prot[k].chain, prot[k].resi, prot[k].icode, prot[k].res, pr[k].q > 0 ? 1 : -1}].push_back(k);
        for (auto& [key, idx] : mq) {
            QGroup g; g.q = get<4>(key); g.idx = idx;
            for (int k : idx) for (int i = 0; i < 3; i++) g.c[i] += pos(prot[k])[i] / idx.size();
            pQ.push_back(g);
        }
    }

    // anelli aromatici e cationi della proteina
    vector<Ring> pRings;
    struct Cation { int idx; string pres; };
    vector<Cation> pCations;
    {
        map<tuple<string,int,char>, map<string,int>> ra;
        for (size_t i = 0; i < prot.size(); i++) ra[{prot[i].chain, prot[i].resi, prot[i].icode}][prot[i].name] = i;
        for (auto& [k, am] : ra) {
            const Atom& any = prot[am.begin()->second];
            const string& R = any.res;
            const map<string,int>& atomsOfRes = am;   // le lambda non possono catturare i structured binding (C++17)
            auto tryRing = [&](vector<string> names, const string& suffix) {
                vector<const Atom*> pts;
                for (auto& n : names) { auto it = atomsOfRes.find(n); if (it == atomsOfRes.end()) return; pts.push_back(&prot[it->second]); }
                Ring rg;
                if (makeRing(pts, resLabel(any) + ":" + suffix, resLabel(any), 0.25, rg)) pRings.push_back(rg);
            };
            if (R == "PHE" || R == "TYR") tryRing({"CG","CD1","CE1","CZ","CE2","CD2"}, "anello");
            if (R == "HIS") tryRing({"CG","ND1","CE1","NE2","CD2"}, "anello");
            if (R == "TRP") {
                tryRing({"CG","CD1","NE1","CE2","CD2"}, "anello5");
                tryRing({"CD2","CE2","CZ2","CH2","CZ3","CE3"}, "anello6");
            }
            if (isNuc(R)) {
                char b = R.back();
                bool pur = b == 'A' || b == 'G' || b == 'I';
                tryRing({"N1","C2","N3","C4","C5","C6"}, pur ? "base6" : "base");
                if (pur) tryRing({"C4","C5","N7","C8","N9"}, "base5");
            }
            if (R == "LYS" && am.count("NZ")) pCations.push_back({am["NZ"], resLabel(any)});
            if (R == "ARG" && am.count("CZ")) pCations.push_back({am["CZ"], resLabel(any)});
        }
    }

    double protB = 0;
    for (auto& a : prot) protB += a.bf;
    protB /= max<size_t>(1, prot.size());
    double scanR = min(6.0, max({P.hbMax, P.hbLong, P.hydMax, P.saltMax, P.repMax, P.halMax, P.metalMax, 3.6}));
    if (hKnownProt) res.info.push_back("Proteina: " + to_string(nHprot) + " idrogeni espliciti nel PDB, ruoli donatore/accettore raffinati.");
    if (!templates.empty()) res.info.push_back("Template caricati: " + to_string(templates.size()));

    // ligandi composti: residui HETATM legati covalentemente (glicani, dipeptidi) formano un unico ligando, come in PLIP
    struct Work { LKey key; vector<Atom> L; set<LKey> comp; };
    vector<Work> works;
    {
        set<LKey> done;
        auto linked = [&](const vector<Atom>& A, const vector<Atom>& B) {
            for (auto& x : A) {
                if (isMetal(x.el)) continue;
                for (auto& y : B) {
                    if (isMetal(y.el)) continue;
                    if (dist(x, y) < covr(x.el) + covr(y.el) + 0.4) return true;
                }
            }
            return false;
        };
        for (auto& [key0, L0] : ligs) {
            if (done.count(key0)) continue;
            Work w; w.key = key0; w.L = L0; w.comp.insert(key0);
            if (mergeLinked) {
                deque<LKey> q; q.push_back(key0);
                while (!q.empty()) {
                    LKey c = q.front(); q.pop_front();
                    for (auto& [k2, at2] : hetRes) {
                        if (w.comp.count(k2) || !linked(hetRes.at(c), at2)) continue;
                        w.comp.insert(k2); q.push_back(k2);
                        w.L.insert(w.L.end(), at2.begin(), at2.end());
                    }
                }
            }
            for (auto& k : w.comp) done.insert(k);
            works.push_back(w);
        }
    }

    for (auto& wk : works) {
        const LKey& key = wk.key;
        vector<Atom>& L = wk.L;
        vector<char> excl(prot.size(), 0);   // atomi di proteina che in realta' sono residui del ligando composto
        for (size_t k = 0; k < prot.size(); k++)
            if (prot[k].het && wk.comp.count(LKey(prot[k].chain, prot[k].resi, prot[k].icode))) excl[k] = 1;
        // connettivita': dal template (per coordinate) oppure geometrica
        string annot;
        bool tmplOK = false, tmplH = false;
        for (auto& t : templates) {
            vector<Atom> cp = L;
            if (applyTemplate(t, cp, annot, tmplH)) { L = cp; tmplOK = true; break; }
        }
        bool hKnown = tmplOK && tmplH;
        if (!tmplOK) {
            for (size_t i = 0; i < L.size(); i++)
                for (size_t j = i + 1; j < L.size(); j++)
                    if (dist(L[i], L[j]) < covr(L[i].el) + covr(L[j].el) + 0.4) {
                        L[i].nb.push_back(j); L[j].nb.push_back(i);
                    }
            auto ih = ligH.find(key);
            if (ih != ligH.end()) {
                int n = 0;
                for (auto& h : ih->second) {
                    int best = -1; double bd = 1.3;
                    for (size_t i = 0; i < L.size(); i++) { double d = nrm(sub(pos(L[i]), h)); if (d < bd) { bd = d; best = i; } }
                    if (best >= 0) { L[best].hs.push_back(h); n++; }
                }
                if (n > 0) { hKnown = true; annot = "idrogeni del PDB (" + to_string(n) + " H), connettivita' geometrica"; }
            }
            if (!hKnown) annot = templates.empty() ? "euristica geometrica" : "euristica geometrica (nessun template corrisponde)";
        }
        string lid = resLabel(L[0]);
        LigandResult lres;
        if (wk.comp.size() > 1) {
            string names;
            for (auto& k : wk.comp) names += (names.empty() ? "" : " + ") + get<0>(k) + ":" + to_string(get<1>(k));
            lres.info.push_back("Ligando composto da " + to_string(wk.comp.size()) + " residui legati (" + names + ")");
        }
        lres.info.push_back("Annotazione ligando " + lid + ": " + annot);
        {
            double lb = 0;
            for (auto& a : L) lb += a.bf;
            lb /= L.size();
            char bb[120]; snprintf(bb, sizeof bb, "B-factor medio: ligando %.1f, resto della struttura %.1f%s", lb, protB,
                                   lb > 2.0 * protB ? "  (attenzione: ligando poco definito)" : "");
            lres.info.push_back(bb);
        }

        vector<Role> lroles(L.size());
        for (size_t i = 0; i < L.size(); i++) lroles[i] = ligandRole(L, i, hKnown);
        vector<vector<int>> lmet(L.size());  // metalli coordinati da ciascun atomo O/N del ligando
        for (size_t i = 0; i < L.size(); i++) {
            if (!isNO(L[i].el)) continue;
            vector<int> cm;
            grid.near(pos(L[i]), 3.5, cm);
            for (int k : cm)
                if (pr[k].metal && dist(L[i], prot[k]) <= metalLimit(prot[k].el, L[i].el)) lmet[i].push_back(k);
            // O/N coordinati a un metallo sono di norma deprotonati; si rispetta pero' un H esplicito (template/PDB)
            if (!metalDonors && !hKnown && !lmet[i].empty()) lroles[i].don = false;
        }

        // metalli che fanno parte del ligando (es. Fe dell'eme): gli N/O coordinati hanno la coppia solitaria impegnata
        vector<char> pcoord(prot.size(), 0);   // atomi della proteina coordinati a un metallo del ligando
        for (size_t j = 0; j < L.size(); j++) {
            if (!lroles[j].metal) continue;
            for (size_t i = 0; i < L.size(); i++) {
                if (!isNO(L[i].el) || dist(L[i], L[j]) > metalLimit(L[j].el, L[i].el)) continue;
                if (!metalDonors && !hKnown) lroles[i].don = false;
                if (L[i].el == "N") lroles[i].acc = false;
            }
            vector<int> cm;
            grid.near(pos(L[j]), 3.0, cm);
            for (int k : cm)
                if (!excl[k] && isNO(prot[k].el) && dist(L[j], prot[k]) <= metalLimit(L[j].el, prot[k].el)) pcoord[k] = 1;
        }

        vector<Hit> hits;
        char buf[160];
        vector<int> cand;

        // --- interazioni atomo-atomo ---
        for (size_t i = 0; i < L.size(); i++) {
            const Atom& l = L[i];
            const Role& lr = lroles[i];
            grid.near(pos(l), scanR, cand);
            for (int k : cand) {
                if (excl[k]) continue;
                const Atom& p = prot[k];
                double d = dist(l, p);
                Role r = pr[k];
                if (pcoord[k]) { r.don = false; if (p.el == "N") r.acc = false; }
                auto add = [&](const string& t, bool fav, const string& note = "") {
                    hits.push_back({t, label(l), label(p), resLabel(p), note, d, fav, pos(l), pos(p)});
                };
                if (lr.metal || r.metal) {
                    if (lr.metal && r.metal) continue;
                    const Atom& other = lr.metal ? p : l;
                    const Atom& metal = lr.metal ? l : p;
                    const Role& orole = lr.metal ? pr[k] : lr;   // ruoli originali: la coordinazione usa la coppia solitaria
                    bool donor = other.el == "S" || (isNO(other.el) && orole.acc);  // serve una coppia solitaria
                    if (donor && d <= min(P.metalMax, metalLimit(metal.el, other.el)))
                        add("Coordinazione metallica", true, d < 1.9 ? "distanza molto corta (modello a bassa risoluzione?)" : "");
                    continue;
                }
                if (d < 2.0) { add("Possibile legame covalente", true); continue; }

                int qq = lr.q * r.q;
                bool sb = P.saltCentroid <= 0 && qq < 0 && d <= P.saltMax;   // con --salt-centroid i ponti salini si valutano per gruppo
                if (qq > 0 && d <= P.repMax) add("Repulsione elettrostatica", false);

                bool polar = isNO(l.el) && isNO(p.el);
                bool hb = false; string hnote;
                if (polar && qq <= 0 && d >= P.hbMin && d <= max(P.hbMax, P.hbLong)) {
                    bool isLong = d > P.hbMax;
                    if (lr.don && r.acc) {
                        auto dr = isLong ? donorTestLong(l, l.nb, L, p) : donorTest(l, l.nb, L, p);
                        if (dr.ok && angleOk(p, pnb(k), prot, l)) { hb = true; hnote = dr.note; }
                    }
                    if (!hb && lr.acc && r.don) {
                        auto dr = isLong ? donorTestLong(p, pnb(k), prot, l) : donorTest(p, pnb(k), prot, l);
                        if (dr.ok && angleOk(l, l.nb, L, p)) { hb = true; hnote = dr.note; }
                    }
                }
                if (sb && hb) add("Ponte salino", true, "anche legame a idrogeno" + (hnote.empty() ? "" : "; " + hnote));
                else if (sb) add("Ponte salino", true);
                else if (hb) add("Legame a idrogeno", true, qq < 0 ? "assistito da carica" + (hnote.empty() ? "" : "; " + hnote) : hnote);

                if (lr.hal && isNOS(p.el) && d <= P.halMax && !l.nb.empty() &&
                    angleDeg(L[l.nb[0]], l, p) >= P.halAngle) add("Legame ad alogeno", true);

                if (lr.hyd && r.hyd && d >= P.hydMin && d <= P.hydMax) add("Contatto idrofobico", true);
                else if (backboneHyd && lr.hyd && p.el == "C" && (p.name == "CA" || p.name == "C") && !isNuc(p.res) &&
                         d >= P.hydMin && d <= P.hydMax) add("Contatto apolare (backbone)", true);

                double thr = polar ? 2.2 : vdw(l.el) + vdw(p.el) - P.clash;
                if (d < thr) {
                    bool shell = false;  // stessa sfera di coordinazione di un metallo: non e' un vero clash
                    for (int m : lmet[i]) if (dist(prot[m], p) <= 3.5) shell = true;
                    if (!shell) add("Clash sterico", false);
                }
            }
        }

        // --- ponti salini tra centroidi dei gruppi carichi (stessa idea di PLIP) ---
        if (P.saltCentroid > 0) {
            vector<int> par(L.size());
            for (size_t i = 0; i < L.size(); i++) par[i] = i;
            function<int(int)> fnd = [&](int x) { return par[x] == x ? x : par[x] = fnd(par[x]); };
            auto sign = [&](size_t i) { return lroles[i].metal ? 0 : (lroles[i].q > 0) - (lroles[i].q < 0); };
            for (size_t i = 0; i < L.size(); i++)
                for (size_t j = i + 1; j < L.size(); j++) {
                    if (!sign(i) || sign(i) != sign(j)) continue;
                    bool joined = find(L[i].nb.begin(), L[i].nb.end(), (int)j) != L[i].nb.end();
                    for (int a : L[i].nb) if (find(L[j].nb.begin(), L[j].nb.end(), a) != L[j].nb.end()) joined = true;
                    if (joined) par[fnd(i)] = fnd(j);
                }
            map<int, QGroup> lg;
            for (size_t i = 0; i < L.size(); i++) if (sign(i)) { QGroup& g = lg[fnd(i)]; g.q = sign(i); g.idx.push_back(i); }
            for (auto& [root, g] : lg) {
                for (int i : g.idx) for (int c = 0; c < 3; c++) g.c[c] += pos(L[i])[c] / g.idx.size();
                for (auto& pg : pQ) {
                    if (g.q * pg.q >= 0) continue;
                    double dc = nrm(sub(g.c, pg.c));
                    if (dc > P.saltCentroid) continue;
                    int k0 = pg.idx[0];
                    if (excl[k0]) continue;
                    hits.push_back({"Ponte salino", label(L[g.idx[0]]), label(prot[k0]), resLabel(prot[k0]),
                                    "centroidi dei gruppi carichi", dc, true, g.c, pg.c});
                }
            }
        }

        // --- anelli: pi-stacking e cation-pi ---
        vector<Ring> lRings;
        for (auto& path : findRings(L)) {
            vector<const Atom*> pts;
            for (int i : path) pts.push_back(&L[i]);
            Ring rg;
            if (makeRing(pts, lid + ":anello(" + L[path[0]].name + ")", lid, 0.15, rg)) lRings.push_back(rg);
        }
        {
            map<pair<string,string>, Hit> bestStack;  // per (tipo, anello proteina) si tiene l'anello del ligando piu' vicino
            for (auto& lg : lRings) {
                for (auto& pg : pRings) {
                    V3 v = sub(pg.c, lg.c);
                    double d = nrm(v);
                    if (d > P.piMax) continue;
                    double ang = acos(min(1.0, fabs(dot(lg.n, pg.n)))) * 180.0 / M_PI;
                    double hL = fabs(dot(v, lg.n)), hP = fabs(dot(v, pg.n));
                    double offL = sqrt(max(0.0, d*d - hL*hL)), offP = sqrt(max(0.0, d*d - hP*hP));
                    double off = min(offL, offP), h = 0.5 * (hL + hP);
                    string type;
                    if (ang <= P.piAngle && h >= 3.0 && h <= 4.5 && off <= P.piOffset) type = "pi-stacking (parallelo)";
                    else if (ang >= 60.0 && off <= 2.0) type = "pi-stacking (a T)";
                    else continue;
                    snprintf(buf, sizeof buf, "angolo normali %.0f; distanza piani %.1f; offset %.1f A", ang, h, off);
                    Hit h2{type, lg.label, pg.label, pg.pres, buf, d, true, lg.c, pg.c};
                    auto k = make_pair(type, pg.label);
                    auto it = bestStack.find(k);
                    if (it == bestStack.end() || d < it->second.d) bestStack[k] = h2;
                }
            }
            for (auto& kv : bestStack) hits.push_back(kv.second);
        }
        auto catpi = [&](const Atom& cat, const Ring& rg, const string& la, const string& pa_, const string& pres, bool ringIsLigand) {
            V3 v = sub(pos(cat), rg.c);
            double d = nrm(v);
            if (d > P.cpMax || d < 1.0) return;
            double ang = acos(min(1.0, fabs(dot(v, rg.n)) / d)) * 180.0 / M_PI;
            if (ang > 30.0) return;
            snprintf(buf, sizeof buf, "angolo asse-normale %.0f", ang);
            hits.push_back({"Cation-pi", la, pa_, pres, buf, d, true, ringIsLigand ? rg.c : pos(cat), ringIsLigand ? pos(cat) : rg.c});
        };
        for (auto& lg : lRings)
            for (auto& c : pCations) catpi(prot[c.idx], lg, lg.label, label(prot[c.idx]), c.pres, true);
        for (size_t i = 0; i < L.size(); i++) {
            if (lroles[i].q <= 0 || lroles[i].metal) continue;
            for (auto& pg : pRings) catpi(L[i], pg, label(L[i]), pg.label, pg.pres, false);
        }

        // --- ponti mediati da acqua ---
        for (size_t i = 0; i < L.size(); i++) {
            const Atom& l = L[i];
            if (!isNO(l.el) || !(lroles[i].don || lroles[i].acc)) continue;
            for (auto& w : waters) {
                double dw = dist(l, w);
                if (dw < P.hbMin || dw > P.hbMax) continue;
                grid.near(pos(w), P.hbMax, cand);
                for (int k : cand) {
                    if (excl[k]) continue;
                    const Atom& p = prot[k];
                    double dp = dist(w, p);
                    if (dp < P.hbMin || !isNO(p.el) || !(pr[k].don || pr[k].acc)) continue;
                    if (dist(l, p) <= P.hbMax) continue;  // gia' contatto diretto
                    double wa = angleDeg(l, w, p);
                    if (wa < P.wbMin || wa > P.wbMax) continue;   // acqua non ponte realistico (quasi collineare o troppo chiusa)
                    snprintf(buf, sizeof buf, "via %s; d(H2O-prot) %.2f", label(w).c_str(), dp);
                    hits.push_back({"Ponte d'acqua", label(l), label(p), resLabel(p), buf, dw, true, pos(l), pos(p), pos(w), true});
                }
            }
        }


        sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
            return tie(a.type, a.la, a.d) < tie(b.type, b.la, b.d);
        });

        // ---- risultato per questo ligando
        lres.id = lid;
        lres.resName = L[0].res;
        lres.chain = L[0].chain;
        lres.resi = L[0].resi;
        lres.annotation = annot;
        lres.nRings = lRings.size();
        if (wk.comp.size() > 1)
            for (auto& k : wk.comp) lres.residues.push_back(get<0>(k) + ":" + to_string(get<1>(k)));
        for (auto& a : L) lres.atoms.push_back(toRecord(a));
        for (size_t i = 0; i < L.size(); i++)
            for (size_t m = 0; m < L[i].nb.size(); m++)
                if (L[i].nb[m] > (int)i) lres.bonds.push_back({(int)i, L[i].nb[m], L[i].bo.empty() ? 0 : L[i].bo[m]});
        for (auto& h : hits) {
            HitRecord hr;
            hr.type = h.type; hr.ligandAtom = h.la; hr.proteinAtom = h.pa; hr.residue = h.pres; hr.note = h.note;
            hr.distance = h.d; hr.favorable = h.fav; hr.ligandPos = h.p1; hr.proteinPos = h.p2; hr.waterPos = h.pw; hr.hasWater = h.hasW;
            if (h.hasW) lres.waters.push_back(h.pw);
            lres.hits.push_back(hr);
        }
        {   // atomi dei residui di proteina che interagiscono (per la grafica)
            set<string> need;
            for (auto& h : hits) need.insert(h.pres);
            vector<int> sel;
            for (size_t k = 0; k < prot.size(); k++)
                if (!excl[k] && need.count(resLabel(prot[k]))) sel.push_back(k);
            for (int k : sel) lres.contactAtoms.push_back(toRecord(prot[k]));
            for (size_t a = 0; a < sel.size(); a++)
                for (size_t b = a + 1; b < sel.size(); b++)
                    if (dist(prot[sel[a]], prot[sel[b]]) < covr(prot[sel[a]].el) + covr(prot[sel[b]].el) + 0.4)
                        lres.contactBonds.push_back({(int)a, (int)b, 0});
        }

        {
            auto isDonorAtom = [&](const Atom& a, const Role& r) { return a.el == "S" || (isNO(a.el) && r.acc); };
            string out;
            vector<int> c2;
            for (size_t k = 0; k < prot.size(); k++) {
                if (!pr[k].metal) continue;
                vector<pair<double,string>> items;
                bool touch = false;
                for (size_t i = 0; i < L.size(); i++) {
                    double d = dist(L[i], prot[k]);
                    if (isDonorAtom(L[i], lroles[i]) && d <= metalLimit(prot[k].el, L[i].el)) {
                        items.push_back({d, label(L[i])}); touch = true;
                    }
                }
                if (!touch) continue;
                grid.near(pos(prot[k]), 3.0, c2);
                for (int j : c2)
                    if ((int)k != j && !excl[j] && !pr[j].metal && isDonorAtom(prot[j], pr[j]) &&
                        dist(prot[j], prot[k]) <= metalLimit(prot[k].el, prot[j].el))
                        items.push_back({dist(prot[j], prot[k]), label(prot[j])});
                for (auto& w : waters)
                    if (dist(w, prot[k]) <= metalLimit(prot[k].el, "O")) items.push_back({dist(w, prot[k]), label(w)});
                sort(items.begin(), items.end());
                char hd[80]; snprintf(hd, sizeof hd, "  %-10s %zu atomi:", label(prot[k]).c_str(), items.size());
                out += hd;
                for (auto& it : items) { char b2[80]; snprintf(b2, sizeof b2, " %s (%.2f)", it.second.c_str(), it.first); out += b2; }
                out += "\n";
            }
            lres.metalSphere = out;
        }
        res.ligands.push_back(std::move(lres));
    }
    return res;
}

}  // namespace ligint
