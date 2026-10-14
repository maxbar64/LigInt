// Lettura dei template SDF / MOL2 e abbinamento al ligando del PDB (per coordinate).
#include "internal.hpp"

namespace ligint {
using namespace detail;
using namespace std;

bool loadSdf(const string& path, vector<Tmpl>& out, vector<string>& warnings) {
    ifstream f(path);
    if (!f) return false;
    string line;
    while (true) {
        Tmpl t; string l1, l2, l3, cnt;
        if (!getline(f, l1)) break;
        t.name = trim(l1);
        if (!getline(f, l2) || !getline(f, l3) || !getline(f, cnt)) break;
        if (cnt.find("V3000") != string::npos) {
            warnings.push_back("SDF V3000 non supportato (molecola '" + t.name + "' saltata)");
            while (getline(f, line) && line.compare(0, 4, "$$$$") != 0) {}
            continue;
        }
        if (cnt.size() < 6) { while (getline(f, line) && line.compare(0, 4, "$$$$") != 0) {} continue; }
        int na = atoi(cnt.substr(0, 3).c_str()), nb = atoi(cnt.substr(3, 3).c_str());
        bool bad = na <= 0 || na > 100000 || nb < 0 || nb > 200000;
        for (int i = 0; i < na && !bad; i++) {
            if (!getline(f, line) || line.size() < 34) { bad = true; break; }
            TAtom a;
            a.p = {atof(line.substr(0, 10).c_str()), atof(line.substr(10, 10).c_str()), atof(line.substr(20, 10).c_str())};
            if (!isfinite(a.p[0]) || !isfinite(a.p[1]) || !isfinite(a.p[2])) { bad = true; break; }
            a.el = normEl(line.substr(31, 3));
            int cc = line.size() >= 39 ? atoi(line.substr(36, 3).c_str()) : 0;
            static const int cmap[8] = {0, 3, 2, 1, 0, -1, -2, -3};
            if (cc >= 1 && cc <= 7) a.fc = cmap[cc];
            t.at.push_back(a);
        }
        for (int i = 0; i < nb && !bad; i++) {
            if (!getline(f, line) || line.size() < 9) { bad = true; break; }
            int ba = atoi(line.substr(0, 3).c_str()) - 1, bb = atoi(line.substr(3, 3).c_str()) - 1;
            if (ba < 0 || bb < 0 || ba >= (int)t.at.size() || bb >= (int)t.at.size() || ba == bb) { bad = true; break; }
            t.bd.push_back({ba, bb, atoi(line.substr(6, 3).c_str())});
        }
        if (bad) {
            warnings.push_back("SDF: molecola '" + t.name + "' malformata, ignorata");
            while (getline(f, line) && line.compare(0, 4, "$$$$") != 0) {}
            continue;
        }
        while (getline(f, line)) {
            if (line.compare(0, 6, "M  END") == 0) break;
            if (line.compare(0, 6, "M  CHG") == 0) {
                istringstream ss(line.substr(6));
                int n = 0; ss >> n;
                for (int k = 0; k < n; k++) {
                    int idx, v;
                    if (!(ss >> idx >> v)) break;
                    if (idx >= 1 && idx <= (int)t.at.size()) t.at[idx - 1].fc = v;
                }
            }
        }
        while (getline(f, line)) if (line.compare(0, 4, "$$$$") == 0) break;
        if (!t.at.empty()) out.push_back(t);
    }
    return true;
}

bool loadMol2(const string& path, vector<Tmpl>& out, vector<string>& warnings) {
    (void)warnings;
    ifstream f(path);
    if (!f) return false;
    string line, sec;
    Tmpl cur; bool have = false, nextName = false;
    map<string,int> idmap;
    auto flush = [&]() {
        if (have && !cur.at.empty()) {
            // cariche formali dedotte dai tipi SYBYL
            for (auto& a : cur.at) {
                if (a.type == "N.4") a.fc = +1;
                if (a.type == "O.co2") a.fc = -1;
            }
            vector<vector<int>> adj(cur.at.size());
            for (auto& b : cur.bd) { adj[b.a].push_back(b.b); adj[b.b].push_back(b.a); }
            for (size_t i = 0; i < cur.at.size(); i++)
                if (cur.at[i].type == "C.cat")
                    for (int j : adj[i]) if (cur.at[j].el == "N") cur.at[j].fc = +1;
            out.push_back(cur);
        }
        cur = Tmpl(); have = false; idmap.clear();
    };
    while (getline(f, line)) {
        if (line.compare(0, 9, "@<TRIPOS>") == 0) {
            sec = trim(line.substr(9));
            if (sec == "MOLECULE") { flush(); have = true; nextName = true; }
            continue;
        }
        if (trim(line).empty() || line[0] == '#') continue;
        if (sec == "MOLECULE" && nextName) { cur.name = trim(line); nextName = false; }
        else if (sec == "ATOM") {
            istringstream ss(line);
            string id, name, type; double x, y, z;
            if (!(ss >> id >> name >> x >> y >> z >> type)) continue;
            TAtom a; a.type = type; a.p = {x, y, z};
            a.el = normEl(type.substr(0, type.find('.')));
            idmap[id] = cur.at.size();
            cur.at.push_back(a);
        } else if (sec == "BOND") {
            istringstream ss(line);
            string id, a, b, ty;
            if (!(ss >> id >> a >> b >> ty)) continue;
            if (!idmap.count(a) || !idmap.count(b)) continue;
            int o = ty == "2" ? 2 : ty == "3" ? 3 : ty == "ar" ? 4 : 1;
            cur.bd.push_back({idmap[a], idmap[b], o});
        }
    }
    flush();
    return true;
}


namespace detail {

bool applyTemplate(const Tmpl& t, vector<Atom>& L, string& info, bool& hasH) {
    vector<int> th;
    for (size_t i = 0; i < t.at.size(); i++) if (!isHel(t.at[i].el)) th.push_back(i);
    if (th.size() != L.size()) return false;
    vector<int> t2l(t.at.size(), -1);
    vector<char> used(t.at.size(), 0);
    for (size_t i = 0; i < L.size(); i++) {
        double best = 0.6; int bj = -1;
        for (int j : th) {
            if (used[j] || t.at[j].el != L[i].el) continue;
            double d = nrm(sub(pos(L[i]), t.at[j].p));
            if (d < best) { best = d; bj = j; }
        }
        if (bj < 0) return false;
        used[bj] = 1; t2l[bj] = i;
        L[i].fc = t.at[bj].fc; L[i].ttype = t.at[bj].type;
        L[i].nb.clear(); L[i].bo.clear(); L[i].hs.clear();
    }
    int nH = 0;
    for (auto& b : t.bd) {
        bool ha = isHel(t.at[b.a].el), hb = isHel(t.at[b.b].el);
        if (!ha && !hb) {
            int p = t2l[b.a], q = t2l[b.b];
            if (p < 0 || q < 0) continue;
            L[p].nb.push_back(q); L[p].bo.push_back(b.order);
            L[q].nb.push_back(p); L[q].bo.push_back(b.order);
        } else if (ha != hb) {
            int hv = t2l[ha ? b.b : b.a];
            if (hv >= 0) { L[hv].hs.push_back(t.at[ha ? b.a : b.b].p); nH++; }
        }
    }
    hasH = nH > 0;
    info = "template '" + t.name + "': " + to_string(th.size()) + " atomi pesanti, " + to_string(nH) + " H" +
           (hasH ? "" : " (senza H: donatori/accettori da euristica)");
    return true;
}

}  // namespace detail
}  // namespace ligint
