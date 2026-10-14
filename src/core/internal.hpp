// Internal details of the library (not part of the public API).
#pragma once

#include "ligint/ligint.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ligint {
namespace detail {
using namespace std;
using ligint::V3;
using ligint::Params;
using ligint::Tmpl;
using ligint::TAtom;
using ligint::TBond;

struct Atom {
    string name, res, chain, el, ttype;
    char icode = ' ';
    int resi = 0, fc = 0;
    bool het = false;   // HETATM record
    double x = 0, y = 0, z = 0, bf = 0;
    vector<int> nb;   // bonds (ligand only)
    vector<int> bo;   // bond order, parallel to nb (empty = unknown; 4 = aromatic)
    vector<V3> hs;    // bonded explicit hydrogens
};
struct Role { bool don = false, acc = false, hyd = false, metal = false, hal = false; int q = 0; };
// Interaction found (internal version of HitRecord, with the coordinates for the graphics).
struct Hit {
    string type, la, pa, pres, note;
    double d;
    bool fav;
    V3 p1{0, 0, 0}, p2{0, 0, 0}, pw{0, 0, 0};
    bool hasW = false;
    int rq = 0;       // charge of the protein atom (salt bridges, repulsions, charge-assisted H-bonds)
    char hbr = 0;     // H-bonds: 'D' protein donates, 'A' protein accepts
};
struct Ring { V3 c, n; string label, pres; };

inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
    return a == string::npos ? "" : s.substr(a, b - a + 1);
}
inline V3 pos(const Atom& a) { return {a.x, a.y, a.z}; }
inline V3 sub(const V3& a, const V3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
inline double dot(const V3& a, const V3& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline double nrm(const V3& a) { return sqrt(dot(a, a)); }
inline double dist(const Atom& a, const Atom& b) { return nrm(sub(pos(a), pos(b))); }
inline double angV(const V3& a, const V3& v, const V3& c) {  // angle a-v-c in degrees
    V3 u = sub(a, v), w = sub(c, v);
    double den = nrm(u) * nrm(w);
    if (den < 1e-9) return 180.0;
    return acos(max(-1.0, min(1.0, dot(u, w) / den))) * 180.0 / M_PI;
}
inline double angleDeg(const Atom& a, const Atom& v, const Atom& c) { return angV(pos(a), pos(v), pos(c)); }
inline bool isMetal(const string& e) {
    static const set<string> m = {"Zn","Mg","Ca","Fe","Mn","Cu","Ni","Co","Cd","Hg","Na","K"};
    return m.count(e);
}
inline double metalLimit(const string& m, const string& o) {
    if (m == "Mg") return 2.5;
    if (m == "Zn") return o == "S" ? 2.7 : 2.6;
    if (m == "Ca") return 2.8;
    if (m == "Na") return 3.0;
    if (m == "K") return 3.4;
    return 2.7;
}
inline bool isNO(const string& e)  { return e == "N" || e == "O"; }
inline bool isNOS(const string& e) { return e == "N" || e == "O" || e == "S"; }
inline double vdw(const string& e) {
    static const map<string,double> v = {{"C",1.70},{"N",1.55},{"O",1.52},{"S",1.80},{"P",1.80},
        {"F",1.47},{"Cl",1.75},{"Br",1.85},{"I",1.98}};
    auto it = v.find(e); return it == v.end() ? 1.8 : it->second;
}
inline double covr(const string& e) {
    static const map<string,double> v = {{"C",.76},{"N",.71},{"O",.66},{"S",1.05},{"P",1.07},
        {"F",.57},{"Cl",1.02},{"Br",1.20},{"I",1.39}};
    auto it = v.find(e); return it == v.end() ? 0.8 : it->second;
}
inline string normEl(string e) {
    e = trim(e);
    if (e.empty()) return e;
    e[0] = toupper(e[0]);
    for (size_t i = 1; i < e.size(); i++) e[i] = tolower(e[i]);
    return e;
}
inline string resLabel(const Atom& a) {
    return a.res + " " + to_string(a.resi) + (a.icode != ' ' ? string(1, a.icode) : "") + a.chain;
}
inline string label(const Atom& a) { return resLabel(a) + ":" + a.name; }
inline bool isHel(const string& e) { return e == "H" || e == "D"; }

// ---------- spatial grid ----------
struct Grid {
    static constexpr double CS = Params::kMaxGridRadius;
    unordered_map<long long, vector<int>> cells;
    const vector<Atom>* A = nullptr;
    static long long key(int i, int j, int k) {
        return ((long long)(i + 1024) << 42) | ((long long)(j + 1024) << 21) | (long long)(k + 1024);
    }
    void build(const vector<Atom>& a) {
        A = &a; cells.clear();
        for (size_t i = 0; i < a.size(); i++)
            cells[key((int)floor(a[i].x / CS), (int)floor(a[i].y / CS), (int)floor(a[i].z / CS))].push_back(i);
    }
    void near(const V3& p, double r, vector<int>& out) const {  // r <= CS
        out.clear();
        int ci = (int)floor(p[0] / CS), cj = (int)floor(p[1] / CS), ck = (int)floor(p[2] / CS);
        for (int i = ci - 1; i <= ci + 1; i++) for (int j = cj - 1; j <= cj + 1; j++)
            for (int k = ck - 1; k <= ck + 1; k++) {
                auto it = cells.find(key(i, j, k));
                if (it == cells.end()) continue;
                for (int idx : it->second) if (nrm(sub(pos((*A)[idx]), p)) <= r) out.push_back(idx);
            }
    }
};

// ---------- functions shared between source files ----------
// templates.cpp
bool applyTemplate(const Tmpl& t, vector<Atom>& L, string& info, bool& hasH);
// chemistry.cpp
bool makeRing(const vector<const Atom*>& pts, const string& lab, const string& pres, double tol, Ring& out);
vector<vector<int>> findRings(const vector<Atom>& L);
bool isNuc(const string& R);
Role proteinRole(const Atom& a);
void refineProteinRoleWithH(const Atom& a, Role& r);
Role ligandRole(const vector<Atom>& L, int i, bool hKnown);

}  // namespace detail
}  // namespace ligint
