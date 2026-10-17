// Cartoon: a smooth ribbon along the C-alpha trace. Loops are round tubes, helices and strands flat ribbons (strands end
// with an arrow). The ribbon is oriented with the C=O direction of the peptide plane, as in the usual cartoon pictures.
#include <algorithm>
#include <cmath>

#include "geom.hpp"
#include "ligint/scene3d.hpp"

namespace ligint {
using namespace geom;

void Mesh::append(const Mesh& o) {
    const int base = static_cast<int>(pos.size());
    pos.insert(pos.end(), o.pos.begin(), o.pos.end());
    nrm.insert(nrm.end(), o.nrm.begin(), o.nrm.end());
    tag.insert(tag.end(), o.tag.begin(), o.tag.end());
    tri.reserve(tri.size() + o.tri.size());
    for (const Tri& t : o.tri) tri.push_back({t[0] + base, t[1] + base, t[2] + base});
}

double meshArea(const Mesh& m) {
    double a = 0;
    for (const Tri& t : m.tri) {
        const V3 p0{m.pos[t[0]][0], m.pos[t[0]][1], m.pos[t[0]][2]};
        const V3 p1{m.pos[t[1]][0], m.pos[t[1]][1], m.pos[t[1]][2]};
        const V3 p2{m.pos[t[2]][0], m.pos[t[2]][1], m.pos[t[2]][2]};
        a += 0.5 * len(cross(sub(p1, p0), sub(p2, p0)));
    }
    return a;
}

Mesh clipMeshNear(const Mesh& m, const std::vector<V3>& points, double radius) {
    Mesh out;
    std::vector<int> remap(m.pos.size(), -1);
    const double r2 = radius * radius;
    for (const Tri& t : m.tri) {
        V3 c{0, 0, 0};
        for (int k = 0; k < 3; k++) for (int i = 0; i < 3; i++) c[i] += m.pos[t[k]][i] / 3.0;
        bool keep = false;
        for (const V3& p : points) {
            const V3 d = sub(c, p);
            if (dot(d, d) <= r2) { keep = true; break; }
        }
        if (!keep) continue;
        Tri nt;
        for (int k = 0; k < 3; k++) {
            int& id = remap[t[k]];
            if (id < 0) {
                id = static_cast<int>(out.pos.size());
                out.pos.push_back(m.pos[t[k]]);
                out.nrm.push_back(m.nrm[t[k]]);
                out.tag.push_back(m.tag[t[k]]);
            }
            nt[k] = id;
        }
        out.tri.push_back(nt);
    }
    return out;
}

namespace {

struct Frame { V3 p, t, w, n; };   // position, tangent, width direction, thickness direction

double smoothstep(double x) { x = std::clamp(x, 0.0, 1.0); return x * x * (3 - 2 * x); }

V3 catmullRom(const V3& p0, const V3& p1, const V3& p2, const V3& p3, double s) {
    const double s2 = s * s, s3 = s2 * s;
    V3 r{0, 0, 0};
    for (int i = 0; i < 3; i++)
        r[i] = 0.5 * ((2 * p1[i]) + (-p0[i] + p2[i]) * s + (2 * p0[i] - 5 * p1[i] + 4 * p2[i] - p3[i]) * s2 + (-p0[i] + 3 * p1[i] - 3 * p2[i] + p3[i]) * s3);
    return r;
}
V3 catmullRomTangent(const V3& p0, const V3& p1, const V3& p2, const V3& p3, double s) {
    const double s2 = s * s;
    V3 r{0, 0, 0};
    for (int i = 0; i < 3; i++)
        r[i] = 0.5 * ((-p0[i] + p2[i]) + 2 * (2 * p0[i] - 5 * p1[i] + 4 * p2[i] - p3[i]) * s + 3 * (-p0[i] + 3 * p1[i] - 3 * p2[i] + p3[i]) * s2);
    return r;
}

struct Profile { double a, b, e; };   // half width, half thickness, exponent of the superellipse

}  // namespace

Mesh buildCartoon(const std::vector<BackboneResidue>& res, const CartoonOptions& opt) {
    Mesh mesh;
    const int N = static_cast<int>(res.size());
    if (N < 2) return mesh;

    // level of detail: the number of triangles stays reasonable for very large proteins
    int subdiv = 8, sides = 16;
    if (opt.quality <= 0) { subdiv = 3; sides = 8; }
    else if (opt.quality == 1) { subdiv = 5; sides = 12; }
    if (N > 2500) { subdiv = std::min(subdiv, 5); sides = std::min(sides, 12); }
    if (N > 8000) { subdiv = std::min(subdiv, 3); sides = std::min(sides, 8); }

    const Profile coil{opt.tubeOnly ? opt.tubeRadius : opt.coilRadius, opt.tubeOnly ? opt.tubeRadius : opt.coilRadius, 2.0};
    const Profile helix{opt.helixHalfWidth, opt.helixHalfThickness, 2.4};
    const Profile strand{opt.strandHalfWidth, opt.strandHalfThickness, 3.2};
    auto profileOf = [&](char ss) {
        if (opt.tubeOnly) return coil;
        return ss == 'H' ? helix : (ss == 'E' ? strand : coil);
    };

    int start = 0;
    while (start < N) {
        int end = start + 1;
        while (end < N && res[end].linkedToPrevious) end++;
        const int n = end - start;             // residues of this segment: [start, end)
        if (n >= 2) {
            // ---- control points and guide vectors
            std::vector<V3> P(n), T(n), G(n);
            for (int k = 0; k < n; k++) P[k] = res[start + k].ca;
            for (int k = 0; k < n; k++) {
                const V3 a = P[std::max(0, k - 1)], b = P[std::min(n - 1, k + 1)];
                T[k] = normalized(sub(b, a));
            }
            V3 prev{0, 0, 0};
            for (int k = 0; k < n; k++) {
                const BackboneResidue& r = res[start + k];
                V3 g = (r.hasO && r.hasC) ? sub(r.o, r.c) : V3{0, 0, 0};
                g = sub(g, mul(T[k], dot(g, T[k])));
                if (len(g) < 0.2) g = k > 0 ? prev : perpendicular(T[k]);
                g = normalized(g);
                if (k > 0 && dot(g, prev) < 0) g = mul(g, -1.0);     // the C=O of a strand alternates: keep one side
                G[k] = prev = g;
            }
            for (int pass = 0; pass < 2; pass++) {                    // smooth the guides (helices and strands)
                std::vector<V3> g2 = G;
                for (int k = 0; k < n; k++) {
                    const char ss = res[start + k].ss;
                    if (ss == 'C' && !opt.tubeOnly) continue;
                    V3 s = mul(G[k], 2.0);
                    if (k > 0) s = add(s, G[k - 1]);
                    if (k + 1 < n) s = add(s, G[k + 1]);
                    g2[k] = normalized(s, G[k]);
                }
                G = g2;
            }

            // arrow heads: the last residue of each strand
            std::vector<char> arrow(n, 0);
            if (!opt.tubeOnly)
                for (int k = 0; k < n; k++)
                    if (res[start + k].ss == 'E' && (k + 1 >= n || res[start + k + 1].ss != 'E') && k >= 1 && res[start + k - 1].ss == 'E') arrow[k] = 1;

            auto profileAt = [&](double u, int k, double s) {          // u = k + s: position along the segment
                // inside an arrow head
                for (int d = -1; d <= 1; d++) {
                    const int L = k + d;
                    if (L < 0 || L >= n || !arrow[L]) continue;
                    const double x = u - (L - 0.6);                      // 0 at the base of the head, 1 at the tip
                    if (x >= 0 && x < 1.0) {
                        return Profile{opt.arrowHalfWidth * (1 - x) + 0.08 * x, strand.b, strand.e};
                    }
                    if (u >= L + 0.4 && u < L + 1.0 && k == L) {         // after the tip: back to the following shape
                        const Profile next = profileOf(k + 1 < n ? res[start + k + 1].ss : 'C');
                        const double t = smoothstep((u - (L + 0.4)) / 0.6);
                        const Profile from{0.08, strand.b, 2.0};
                        return Profile{from.a + (next.a - from.a) * t, from.b + (next.b - from.b) * t, from.e + (next.e - from.e) * t};
                    }
                }
                const Profile p0 = profileOf(res[start + k].ss);
                const Profile p1 = profileOf(res[start + std::min(n - 1, k + 1)].ss);
                const double t = smoothstep(s);
                return Profile{p0.a + (p1.a - p0.a) * t, p0.b + (p1.b - p0.b) * t, p0.e + (p1.e - p0.e) * t};
            };

            // ---- rings along the spline
            const int segs = n - 1;
            const int rings = segs * subdiv + 1;
            const int base = static_cast<int>(mesh.pos.size());
            V3 prevW{0, 0, 0};
            for (int ri = 0; ri < rings; ri++) {
                const int k = std::min(segs - 1, ri / subdiv);
                const double s = ri == rings - 1 ? 1.0 : static_cast<double>(ri % subdiv) / subdiv;
                const V3 p1 = P[k], p2 = P[k + 1];
                // phantom end points keep the tangent smooth at the ends of the segment
                const V3 c0 = k == 0 ? add(mul(P[0], 2.0), mul(P[1], -1.0)) : P[k - 1];
                const V3 c3 = k + 2 >= n ? add(mul(P[n - 1], 2.0), mul(P[n - 2], -1.0)) : P[k + 2];
                const V3 pos = catmullRom(c0, p1, p2, c3, s);
                const V3 tan = normalized(catmullRomTangent(c0, p1, p2, c3, s), T[k]);
                const double ts = smoothstep(s);
                V3 w = normalized(lerp(G[k], G[k + 1], ts), G[k]);
                w = sub(w, mul(tan, dot(w, tan)));
                w = normalized(w, ri > 0 ? prevW : perpendicular(tan));
                prevW = w;
                const V3 nrm = normalized(cross(tan, w));
                const double u = k + s;
                const Profile pf = profileAt(u, k, s);
                const double tag = start + u;
                for (int j = 0; j < sides; j++) {
                    const double th = 2.0 * M_PI * j / sides;
                    const double cs = std::cos(th), sn = std::sin(th);
                    const double ex = 2.0 / pf.e;
                    const double x = pf.a * (cs < 0 ? -1 : 1) * std::pow(std::fabs(cs), ex);
                    const double y = pf.b * (sn < 0 ? -1 : 1) * std::pow(std::fabs(sn), ex);
                    // gradient of |x/a|^e + |y/b|^e (outward normal)
                    const double gx = pf.a > 1e-6 ? (x < 0 ? -1 : 1) * std::pow(std::fabs(x) / pf.a, pf.e - 1) / pf.a : 0;
                    const double gy = pf.b > 1e-6 ? (y < 0 ? -1 : 1) * std::pow(std::fabs(y) / pf.b, pf.e - 1) / pf.b : 0;
                    const V3 nl = normalized(add(mul(w, gx), mul(nrm, gy)), w);
                    const V3 pt = add(pos, add(mul(w, x), mul(nrm, y)));
                    mesh.pos.push_back({(float)pt[0], (float)pt[1], (float)pt[2]});
                    mesh.nrm.push_back({(float)nl[0], (float)nl[1], (float)nl[2]});
                    mesh.tag.push_back((float)tag);
                }
            }
            for (int ri = 0; ri + 1 < rings; ri++)
                for (int j = 0; j < sides; j++) {
                    const int j2 = (j + 1) % sides;
                    const int a = base + ri * sides + j, b = base + ri * sides + j2;
                    const int c = base + (ri + 1) * sides + j, d = base + (ri + 1) * sides + j2;
                    mesh.tri.push_back({a, b, c});
                    mesh.tri.push_back({b, d, c});
                }
            // end caps
            for (int cap = 0; cap < 2; cap++) {
                const int ring = cap == 0 ? 0 : rings - 1;
                V3 center{0, 0, 0};
                for (int j = 0; j < sides; j++)
                    for (int i = 0; i < 3; i++) center[i] += mesh.pos[base + ring * sides + j][i] / sides;
                const V3 tn = normalized(cap == 0 ? sub(P[0], P[1]) : sub(P[n - 1], P[n - 2]));
                const int ci = static_cast<int>(mesh.pos.size());
                mesh.pos.push_back({(float)center[0], (float)center[1], (float)center[2]});
                mesh.nrm.push_back({(float)tn[0], (float)tn[1], (float)tn[2]});
                mesh.tag.push_back(mesh.tag[base + ring * sides]);
                // the cap needs its own ring of vertices with the flat normal
                const int rb = static_cast<int>(mesh.pos.size());
                for (int j = 0; j < sides; j++) {
                    mesh.pos.push_back(mesh.pos[base + ring * sides + j]);
                    mesh.nrm.push_back({(float)tn[0], (float)tn[1], (float)tn[2]});
                    mesh.tag.push_back(mesh.tag[base + ring * sides + j]);
                }
                for (int j = 0; j < sides; j++) {
                    const int j2 = (j + 1) % sides;
                    if (cap == 1) mesh.tri.push_back({ci, rb + j, rb + j2});
                    else mesh.tri.push_back({ci, rb + j2, rb + j});
                }
            }
        }
        start = end;
    }
    return mesh;
}

}  // namespace ligint
