// Molecular surfaces on a grid.
//
// 1. A signed distance to the solvent-accessible surface (SAS): for every atom, distance to the centre minus (radius + probe).
// 2. The cavities that the solvent cannot reach are filled (flood fill from the outside of the box).
// 3. SAS: iso-surface 0 of that field. Solvent-excluded surface (SES, "Connolly"): the region reachable by the centre of the
//    probe is dilated by the probe radius: iso-surface (probe - distance to that region), by an exact Euclidean distance
//    transform (Felzenszwalb & Huttenlocher).
// 4. The iso-surface is extracted with marching tetrahedra (no ambiguous cases), smoothed (Taubin) and given smooth normals.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

#include "geom.hpp"
#include "ligint/scene3d.hpp"

namespace ligint {
using namespace geom;

namespace {

struct Grid3 {
    int nx = 0, ny = 0, nz = 0;
    double h = 0.5;
    V3 origin{0, 0, 0};
    std::size_t size() const { return static_cast<std::size_t>(nx) * ny * nz; }
    std::size_t id(int i, int j, int k) const { return (static_cast<std::size_t>(k) * ny + j) * nx + i; }
};

void edt1d(const std::vector<double>& f, std::vector<double>& d, int n, std::vector<int>& v, std::vector<double>& z) {
    int k = 0;
    v[0] = 0;
    z[0] = -1e30;
    z[1] = 1e30;
    for (int q = 1; q < n; q++) {
        double s;
        while (true) {
            const int vk = v[k];
            s = ((f[q] + static_cast<double>(q) * q) - (f[vk] + static_cast<double>(vk) * vk)) / (2.0 * q - 2.0 * vk);
            if (s <= z[k]) k--; else break;
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = 1e30;
    }
    k = 0;
    for (int q = 0; q < n; q++) {
        while (z[k + 1] < q) k++;
        const double dq = q - v[k];
        d[q] = dq * dq + f[v[k]];
    }
}

// Squared distance (in voxels) from every voxel to the nearest voxel with mask != 0.
void distanceTransform(const Grid3& g, const std::vector<char>& mask, std::vector<float>& out) {
    out.assign(g.size(), 0.f);
    const double INF = 1e12;
    for (std::size_t i = 0; i < out.size(); i++) out[i] = mask[i] ? 0.f : static_cast<float>(INF);
    const int m = std::max({g.nx, g.ny, g.nz});
    std::vector<double> f(m), d(m), z(m + 1);
    std::vector<int> v(m);
    for (int k = 0; k < g.nz; k++)           // along x
        for (int j = 0; j < g.ny; j++) {
            for (int i = 0; i < g.nx; i++) f[i] = out[g.id(i, j, k)];
            edt1d(f, d, g.nx, v, z);
            for (int i = 0; i < g.nx; i++) out[g.id(i, j, k)] = static_cast<float>(std::min(d[i], INF));
        }
    for (int k = 0; k < g.nz; k++)           // along y
        for (int i = 0; i < g.nx; i++) {
            for (int j = 0; j < g.ny; j++) f[j] = out[g.id(i, j, k)];
            edt1d(f, d, g.ny, v, z);
            for (int j = 0; j < g.ny; j++) out[g.id(i, j, k)] = static_cast<float>(std::min(d[j], INF));
        }
    for (int j = 0; j < g.ny; j++)           // along z
        for (int i = 0; i < g.nx; i++) {
            for (int k = 0; k < g.nz; k++) f[k] = out[g.id(i, j, k)];
            edt1d(f, d, g.nz, v, z);
            for (int k = 0; k < g.nz; k++) out[g.id(i, j, k)] = static_cast<float>(std::min(d[k], INF));
        }
}

// Marching tetrahedra: the six tetrahedra around the diagonal 0-6 of every cell.
constexpr int kTets[6][4] = {{0, 5, 1, 6}, {0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}};
constexpr int kCorner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// One pass of the binomial filter [1 2 1]/4 along each axis. The distance transform is quantized (its values are sqrt(n) voxels):
// without this the iso-surface can lock onto the nodes of the grid and lose part of its area.
void blur121(const Grid3& g, std::vector<float>& f) {
    std::vector<float> t(f.size());
    auto pass = [&](int di, int dj, int dk) {
        for (int k = 0; k < g.nz; k++)
            for (int j = 0; j < g.ny; j++)
                for (int i = 0; i < g.nx; i++) {
                    const int i0 = std::max(0, i - di), i1 = std::min(g.nx - 1, i + di);
                    const int j0 = std::max(0, j - dj), j1 = std::min(g.ny - 1, j + dj);
                    const int k0 = std::max(0, k - dk), k1 = std::min(g.nz - 1, k + dk);
                    t[g.id(i, j, k)] = 0.25f * f[g.id(i0, j0, k0)] + 0.5f * f[g.id(i, j, k)] + 0.25f * f[g.id(i1, j1, k1)];
                }
        f.swap(t);
    };
    pass(1, 0, 0);
    pass(0, 1, 0);
    pass(0, 0, 1);
}

// The distance transform measures from voxel centres, while the region that the probe can reach extends a little beyond
// them, and the blur above moves the surface slightly inwards: 0.04 voxel keeps the area of spheres of the usual atomic radii
// within a few percent of the exact value (calibrated on spheres and on a protein at several grid spacings).
constexpr double kVoxelBias = 0.04;

Mesh isoSurface(const Grid3& g, const std::vector<float>& s) {
    Mesh mesh;
    std::unordered_map<std::uint64_t, int> edgeVertex;
    edgeVertex.reserve(1 << 16);
    std::vector<V3f> grads;                                 // gradient of s at the vertices (orientation of the triangles)

    auto gradAt = [&](int i, int j, int k) {
        auto at = [&](int a, int b, int c) {
            a = std::clamp(a, 0, g.nx - 1); b = std::clamp(b, 0, g.ny - 1); c = std::clamp(c, 0, g.nz - 1);
            return static_cast<double>(s[g.id(a, b, c)]);
        };
        return V3{at(i + 1, j, k) - at(i - 1, j, k), at(i, j + 1, k) - at(i, j - 1, k), at(i, j, k + 1) - at(i, j, k - 1)};
    };
    auto cellOf = [&](std::size_t id, int& i, int& j, int& k) {
        i = static_cast<int>(id % g.nx);
        j = static_cast<int>((id / g.nx) % g.ny);
        k = static_cast<int>(id / (static_cast<std::size_t>(g.nx) * g.ny));
    };
    auto vertexOn = [&](std::size_t a, std::size_t b) {
        if (a > b) std::swap(a, b);
        const std::uint64_t key = (static_cast<std::uint64_t>(a) << 32) ^ static_cast<std::uint64_t>(b);
        auto it = edgeVertex.find(key);
        if (it != edgeVertex.end()) return it->second;
        int ai, aj, ak, bi, bj, bk;
        cellOf(a, ai, aj, ak);
        cellOf(b, bi, bj, bk);
        const double sa = s[a], sb = s[b];
        const double t = std::clamp(sa / (sa - sb), 0.0, 1.0);
        const V3 pa{g.origin[0] + ai * g.h, g.origin[1] + aj * g.h, g.origin[2] + ak * g.h};
        const V3 pb{g.origin[0] + bi * g.h, g.origin[1] + bj * g.h, g.origin[2] + bk * g.h};
        const V3 p = lerp(pa, pb, t);
        const V3 gr = lerp(gradAt(ai, aj, ak), gradAt(bi, bj, bk), t);
        const int idx = static_cast<int>(mesh.pos.size());
        mesh.pos.push_back({(float)p[0], (float)p[1], (float)p[2]});
        grads.push_back({(float)gr[0], (float)gr[1], (float)gr[2]});
        edgeVertex.emplace(key, idx);
        return idx;
    };
    auto emit = [&](int a, int b, int c) {
        if (a == b || b == c || a == c) return;
        const V3 p0{mesh.pos[a][0], mesh.pos[a][1], mesh.pos[a][2]}, p1{mesh.pos[b][0], mesh.pos[b][1], mesh.pos[b][2]},
            p2{mesh.pos[c][0], mesh.pos[c][1], mesh.pos[c][2]};
        const V3 n = cross(sub(p1, p0), sub(p2, p0));
        if (len(n) < 1e-12) return;
        const V3 ref{grads[a][0] + grads[b][0] + grads[c][0], grads[a][1] + grads[b][1] + grads[c][1], grads[a][2] + grads[b][2] + grads[c][2]};
        if (dot(n, ref) < 0) std::swap(b, c);               // s grows outwards: the normals point to the solvent
        mesh.tri.push_back({a, b, c});
    };

    for (int k = 0; k + 1 < g.nz; k++)
        for (int j = 0; j + 1 < g.ny; j++)
            for (int i = 0; i + 1 < g.nx; i++) {
                std::size_t cid[8];
                float cv[8];
                int inside = 0;
                for (int c = 0; c < 8; c++) {
                    cid[c] = g.id(i + kCorner[c][0], j + kCorner[c][1], k + kCorner[c][2]);
                    cv[c] = s[cid[c]];
                    inside += cv[c] < 0.f;
                }
                if (inside == 0 || inside == 8) continue;
                for (const auto& t : kTets) {
                    int in[4], out[4], ni = 0, no = 0;
                    for (int q = 0; q < 4; q++) (cv[t[q]] < 0.f ? in[ni++] : out[no++]) = t[q];
                    if (ni == 0 || no == 0) continue;
                    if (ni == 1) {
                        emit(vertexOn(cid[in[0]], cid[out[0]]), vertexOn(cid[in[0]], cid[out[1]]), vertexOn(cid[in[0]], cid[out[2]]));
                    } else if (no == 1) {
                        emit(vertexOn(cid[out[0]], cid[in[0]]), vertexOn(cid[out[0]], cid[in[1]]), vertexOn(cid[out[0]], cid[in[2]]));
                    } else {
                        const int a = vertexOn(cid[in[0]], cid[out[0]]), b = vertexOn(cid[in[0]], cid[out[1]]);
                        const int c = vertexOn(cid[in[1]], cid[out[1]]), d = vertexOn(cid[in[1]], cid[out[0]]);
                        emit(a, b, c);
                        emit(a, c, d);
                    }
                }
            }
    mesh.nrm.assign(mesh.pos.size(), V3f{0, 0, 1});
    mesh.tag.assign(mesh.pos.size(), 0.f);
    return mesh;
}

void recomputeNormals(Mesh& m) {
    std::vector<V3> acc(m.pos.size(), V3{0, 0, 0});
    for (const Tri& t : m.tri) {
        const V3 p0{m.pos[t[0]][0], m.pos[t[0]][1], m.pos[t[0]][2]}, p1{m.pos[t[1]][0], m.pos[t[1]][1], m.pos[t[1]][2]},
            p2{m.pos[t[2]][0], m.pos[t[2]][1], m.pos[t[2]][2]};
        const V3 n = cross(sub(p1, p0), sub(p2, p0));      // length = 2 * area: area-weighted
        for (int k = 0; k < 3; k++) acc[t[k]] = add(acc[t[k]], n);
    }
    for (std::size_t i = 0; i < acc.size(); i++) {
        const V3 n = normalized(acc[i], V3{0, 0, 1});
        m.nrm[i] = {(float)n[0], (float)n[1], (float)n[2]};
    }
}

// Taubin smoothing (lambda | mu): smooths without shrinking the surface.
void taubinSmooth(Mesh& m, int iterations) {
    if (iterations <= 0 || m.pos.empty()) return;
    const std::size_t n = m.pos.size();
    std::vector<int> start(n + 1, 0);
    for (const Tri& t : m.tri) for (int k = 0; k < 3; k++) start[t[k] + 1] += 2;
    for (std::size_t i = 0; i < n; i++) start[i + 1] += start[i];
    std::vector<int> fill(start.begin(), start.end() - 1), adj(start[n]);
    for (const Tri& t : m.tri)
        for (int k = 0; k < 3; k++) {
            adj[fill[t[k]]++] = t[(k + 1) % 3];
            adj[fill[t[k]]++] = t[(k + 2) % 3];
        }
    std::vector<V3f> tmp(n);
    auto step = [&](float lambda) {
        for (std::size_t i = 0; i < n; i++) {
            const int a = start[i], b = start[i + 1];
            if (a == b) { tmp[i] = m.pos[i]; continue; }
            double c[3] = {0, 0, 0};
            for (int q = a; q < b; q++) for (int d = 0; d < 3; d++) c[d] += m.pos[adj[q]][d];
            for (int d = 0; d < 3; d++) tmp[i][d] = m.pos[i][d] + lambda * static_cast<float>(c[d] / (b - a) - m.pos[i][d]);
        }
        m.pos.swap(tmp);
    };
    for (int it = 0; it < iterations; it++) { step(0.5f); step(-0.53f); }
}

}  // namespace

Mesh buildSurface(const std::vector<ProteinAtom>& atoms, const SurfaceOptions& opt) {
    Mesh empty;
    if (atoms.empty()) return empty;
    const double probe = std::max(0.0, opt.probe);

    // ---- grid
    V3 lo = atoms[0].pos, hi = atoms[0].pos;
    double rmax = 0;
    for (const auto& a : atoms) {
        for (int i = 0; i < 3; i++) { lo[i] = std::min(lo[i], a.pos[i]); hi[i] = std::max(hi[i], a.pos[i]); }
        rmax = std::max<double>(rmax, a.radius);
    }
    Grid3 g;
    double h = std::max(0.2, opt.spacing);
    const double pad = rmax + probe + 3.0 * h;
    for (int attempt = 0; attempt < 8; attempt++) {
        g.h = h;
        g.nx = static_cast<int>(std::ceil((hi[0] - lo[0] + 2 * pad) / h)) + 1;
        g.ny = static_cast<int>(std::ceil((hi[1] - lo[1] + 2 * pad) / h)) + 1;
        g.nz = static_cast<int>(std::ceil((hi[2] - lo[2] + 2 * pad) / h)) + 1;
        if (g.size() <= 7000000ull) break;
        h *= 1.15;
    }
    g.origin = {lo[0] - pad, lo[1] - pad, lo[2] - pad};
    const double BIG = 2.0 * g.h;

    // ---- signed distance to the solvent-accessible surface
    std::vector<float> sas(g.size(), static_cast<float>(BIG));
    for (const auto& a : atoms) {
        const double R = a.radius + probe;
        const double reach = R + BIG;
        const int i0 = std::max(0, static_cast<int>(std::floor((a.pos[0] - reach - g.origin[0]) / g.h)));
        const int i1 = std::min(g.nx - 1, static_cast<int>(std::ceil((a.pos[0] + reach - g.origin[0]) / g.h)));
        const int j0 = std::max(0, static_cast<int>(std::floor((a.pos[1] - reach - g.origin[1]) / g.h)));
        const int j1 = std::min(g.ny - 1, static_cast<int>(std::ceil((a.pos[1] + reach - g.origin[1]) / g.h)));
        const int k0 = std::max(0, static_cast<int>(std::floor((a.pos[2] - reach - g.origin[2]) / g.h)));
        const int k1 = std::min(g.nz - 1, static_cast<int>(std::ceil((a.pos[2] + reach - g.origin[2]) / g.h)));
        for (int k = k0; k <= k1; k++) {
            const double dz = g.origin[2] + k * g.h - a.pos[2];
            for (int j = j0; j <= j1; j++) {
                const double dy = g.origin[1] + j * g.h - a.pos[1];
                const double dyz = dy * dy + dz * dz;
                for (int i = i0; i <= i1; i++) {
                    const double dx = g.origin[0] + i * g.h - a.pos[0];
                    const double d = std::sqrt(dx * dx + dyz) - R;
                    float& v = sas[g.id(i, j, k)];
                    if (d < v) v = static_cast<float>(d);
                }
            }
        }
    }

    // ---- exterior (reachable from the border of the box) and cavities
    std::vector<char> exterior(g.size(), 0);
    {
        std::vector<std::size_t> stack;
        auto push = [&](int i, int j, int k) {
            const std::size_t id = g.id(i, j, k);
            if (!exterior[id] && sas[id] > 0.f) { exterior[id] = 1; stack.push_back(id); }
        };
        for (int k = 0; k < g.nz; k++) for (int j = 0; j < g.ny; j++) { push(0, j, k); push(g.nx - 1, j, k); }
        for (int k = 0; k < g.nz; k++) for (int i = 0; i < g.nx; i++) { push(i, 0, k); push(i, g.ny - 1, k); }
        for (int j = 0; j < g.ny; j++) for (int i = 0; i < g.nx; i++) { push(i, j, 0); push(i, j, g.nz - 1); }
        while (!stack.empty()) {
            const std::size_t id = stack.back();
            stack.pop_back();
            const int i = static_cast<int>(id % g.nx), j = static_cast<int>((id / g.nx) % g.ny), k = static_cast<int>(id / (static_cast<std::size_t>(g.nx) * g.ny));
            if (i > 0) push(i - 1, j, k);
            if (i + 1 < g.nx) push(i + 1, j, k);
            if (j > 0) push(i, j - 1, k);
            if (j + 1 < g.ny) push(i, j + 1, k);
            if (k > 0) push(i, j, k - 1);
            if (k + 1 < g.nz) push(i, j, k + 1);
        }
    }

    std::vector<float> field;
    if (opt.kind == SurfaceKind::SolventAccessible) {
        field = std::move(sas);
        for (std::size_t i = 0; i < field.size(); i++)
            if (!exterior[i] && field[i] > 0.f) field[i] = -0.5f * static_cast<float>(g.h);    // cavity: inside
    } else {
        std::vector<float> d2;
        distanceTransform(g, exterior, d2);
        field.resize(g.size());
        for (std::size_t i = 0; i < field.size(); i++) field[i] = static_cast<float>(probe + kVoxelBias * g.h - std::sqrt(static_cast<double>(d2[i])) * g.h);
        sas.clear();
        sas.shrink_to_fit();
        blur121(g, field);
    }

    Mesh mesh = isoSurface(g, field);
    taubinSmooth(mesh, opt.smoothing);
    recomputeNormals(mesh);

    // ---- nearest atom of every vertex
    {
        const double cell = 5.0;
        std::unordered_map<std::int64_t, std::vector<int>> hash;
        auto ck = [&](int i, int j, int k) { return ((std::int64_t)(i + 8192) << 42) | ((std::int64_t)(j + 8192) << 21) | (std::int64_t)(k + 8192); };
        for (std::size_t i = 0; i < atoms.size(); i++)
            hash[ck((int)std::floor(atoms[i].pos[0] / cell), (int)std::floor(atoms[i].pos[1] / cell), (int)std::floor(atoms[i].pos[2] / cell))].push_back((int)i);
        for (std::size_t v = 0; v < mesh.pos.size(); v++) {
            const V3 p{mesh.pos[v][0], mesh.pos[v][1], mesh.pos[v][2]};
            const int ci = (int)std::floor(p[0] / cell), cj = (int)std::floor(p[1] / cell), ck2 = (int)std::floor(p[2] / cell);
            double best = 1e30;
            int bi = 0;
            for (int a = ci - 1; a <= ci + 1; a++)
                for (int b = cj - 1; b <= cj + 1; b++)
                    for (int c = ck2 - 1; c <= ck2 + 1; c++) {
                        auto it = hash.find(ck(a, b, c));
                        if (it == hash.end()) continue;
                        for (int idx : it->second) {
                            const double d = dist(p, atoms[idx].pos) - atoms[idx].radius;
                            if (d < best) { best = d; bi = idx; }
                        }
                    }
            mesh.tag[v] = static_cast<float>(bi);
        }
    }
    return mesh;
}

}  // namespace ligint
