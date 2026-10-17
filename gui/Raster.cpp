#include "Raster.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>

namespace gui {
namespace raster {
namespace {

constexpr float kFar = std::numeric_limits<float>::max();

template <class F>
void parallelFor(int n, int threads, F&& fn) {
    if (n <= 0) return;
    threads = std::max(1, std::min(threads, n));
    if (threads == 1) {
        for (int i = 0; i < n; i++) fn(i);
        return;
    }
    std::atomic<int> next{0};
    auto worker = [&]() {
        for (;;) {
            const int i = next.fetch_add(1);
            if (i >= n) break;
            fn(i);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(threads - 1);
    for (int t = 1; t < threads; t++) pool.emplace_back(worker);
    worker();
    for (auto& t : pool) t.join();
}

int defaultThreads() {
    const unsigned hw = std::thread::hardware_concurrency();
    return static_cast<int>(std::clamp<unsigned>(hw == 0 ? 2 : hw, 1u, 8u));
}

struct Normalized3 { float x, y, z; };
Normalized3 unit(float x, float y, float z) {
    const float l = std::sqrt(x * x + y * y + z * z);
    return {x / l, y / l, z / l};
}

}  // namespace

// ----------------------------------------------------------------------------------------------
// Camera
// ----------------------------------------------------------------------------------------------
ligint::V3 Camera::toView(const ligint::V3& p) const {
    const double dx = p[0] - center[0], dy = p[1] - center[1], dz = p[2] - center[2];
    return {rot[0] * dx + rot[1] * dy + rot[2] * dz, rot[3] * dx + rot[4] * dy + rot[5] * dz, rot[6] * dx + rot[7] * dy + rot[8] * dz};
}

Camera::Point Camera::project(const ligint::V3& p) const {
    const ligint::V3 v = toView(p);
    const double base = distance > 0 ? distance : 1000.0;
    const double w = base - v[2];
    double f = scale;
    bool ok = true;
    if (distance > 0) {
        ok = w > 0.05 * distance;
        f = scale * distance / std::max(w, 1e-6);
    }
    return {width * 0.5 + panX + f * v[0], height * 0.5 + panY - f * v[1], w, ok};
}

// ----------------------------------------------------------------------------------------------
// Renderer
// ----------------------------------------------------------------------------------------------
Renderer::Renderer(int width, int height, int supersample)
    : width_(std::max(1, width)), height_(std::max(1, height)), ss_(std::max(1, supersample)) {
    // keep the sample buffers within a reasonable size
    while (ss_ > 1 && static_cast<long long>(width_) * height_ * ss_ * ss_ > 24000000LL) ss_--;
    W_ = width_ * ss_;
    H_ = height_ * ss_;
    depth_.assign(static_cast<std::size_t>(W_) * H_, kFar);
    color_.assign(static_cast<std::size_t>(W_) * H_, 0xffffff);
    cam_.width = width_;
    cam_.height = height_;
}

void Renderer::setCamera(const Camera& c) {
    cam_ = c;
    cam_.width = width_;
    cam_.height = height_;
}

void Renderer::setLighting(const Lighting& l) {
    light_ = l;
    std::fill(color_.begin(), color_.end(), l.background);
    std::fill(depth_.begin(), depth_.end(), kFar);
}

int Renderer::bandCount() const { return std::max(1, (threads_ > 0 ? threads_ : defaultThreads()) * 4); }

void Renderer::drawMesh(const MeshView& mv) {
    if (!mv.mesh || !mv.colors || mv.mesh->tri.empty()) return;
    const ligint::Mesh& m = *mv.mesh;
    const std::vector<std::uint32_t>& colors = *mv.colors;
    const int nv = static_cast<int>(m.pos.size());
    const int nt = static_cast<int>(m.tri.size());
    const int nthreads = threads_ > 0 ? threads_ : defaultThreads();

    struct V { float x, y, z, nx, ny, nz, r, g, b; };
    std::vector<V> view(nv);
    const auto& R = cam_.rot;
    parallelFor((nv + 4095) / 4096, nthreads, [&](int chunk) {
        const int a = chunk * 4096, b = std::min(nv, a + 4096);
        for (int i = a; i < b; i++) {
            const double dx = m.pos[i][0] - cam_.center[0], dy = m.pos[i][1] - cam_.center[1], dz = m.pos[i][2] - cam_.center[2];
            V& v = view[i];
            v.x = static_cast<float>(R[0] * dx + R[1] * dy + R[2] * dz);
            v.y = static_cast<float>(R[3] * dx + R[4] * dy + R[5] * dz);
            v.z = static_cast<float>(R[6] * dx + R[7] * dy + R[8] * dz);
            const double nx = m.nrm[i][0], ny = m.nrm[i][1], nz = m.nrm[i][2];
            v.nx = static_cast<float>(R[0] * nx + R[1] * ny + R[2] * nz);
            v.ny = static_cast<float>(R[3] * nx + R[4] * ny + R[5] * nz);
            v.nz = static_cast<float>(R[6] * nx + R[7] * ny + R[8] * nz);
            const std::uint32_t c = colors[i];
            v.r = static_cast<float>((c >> 16) & 255);
            v.g = static_cast<float>((c >> 8) & 255);
            v.b = static_cast<float>(c & 255);
        }
    });

    const bool persp = cam_.distance > 0;
    const double base = persp ? cam_.distance : 1000.0;
    const float zLimit = static_cast<float>(persp ? std::min(cam_.clipZ, base * 0.95) : cam_.clipZ);
    const float sx0 = static_cast<float>(W_ * 0.5 + cam_.panX * ss_), sy0 = static_cast<float>(H_ * 0.5 + cam_.panY * ss_);
    const float scale = static_cast<float>(cam_.scale * ss_);

    const int chunkSize = 2048;
    const int chunks = (nt + chunkSize - 1) / chunkSize;
    std::vector<std::vector<PreTri>> out(chunks);
    parallelFor(chunks, nthreads, [&](int chunk) {
        auto& dst = out[chunk];
        dst.reserve(chunkSize);
        const int a = chunk * chunkSize, b = std::min(nt, a + chunkSize);
        for (int t = a; t < b; t++) {
            V poly[4];
            int np = 0;
            const V* in[3] = {&view[m.tri[t][0]], &view[m.tri[t][1]], &view[m.tri[t][2]]};
            // clip against the plane z = zLimit (keep z <= zLimit)
            for (int k = 0; k < 3; k++) {
                const V& p = *in[k];
                const V& q = *in[(k + 1) % 3];
                const bool pin = p.z <= zLimit, qin = q.z <= zLimit;
                if (pin && np < 4) poly[np++] = p;
                if (pin != qin && np < 4) {
                    const float s = (zLimit - p.z) / (q.z - p.z);
                    V c;
                    c.x = p.x + (q.x - p.x) * s; c.y = p.y + (q.y - p.y) * s; c.z = zLimit;
                    c.nx = p.nx + (q.nx - p.nx) * s; c.ny = p.ny + (q.ny - p.ny) * s; c.nz = p.nz + (q.nz - p.nz) * s;
                    c.r = p.r + (q.r - p.r) * s; c.g = p.g + (q.g - p.g) * s; c.b = p.b + (q.b - p.b) * s;
                    poly[np++] = c;
                }
            }
            if (np < 3) continue;
            Vtx sv[4];
            for (int k = 0; k < np; k++) {
                const V& p = poly[k];
                const float w = static_cast<float>(base - p.z);
                const float f = persp ? static_cast<float>(scale * base / w) : scale;
                sv[k] = {sx0 + f * p.x, sy0 - f * p.y, persp ? 1.0f / w : 1.0f, w, p.nx, p.ny, p.nz, p.r, p.g, p.b};
            }
            for (int k = 1; k + 1 < np; k++) {
                PreTri pt;
                pt.v[0] = sv[0]; pt.v[1] = sv[k]; pt.v[2] = sv[k + 1];
                const float area2 = (pt.v[1].x - pt.v[0].x) * (pt.v[2].y - pt.v[0].y) - (pt.v[1].y - pt.v[0].y) * (pt.v[2].x - pt.v[0].x);
                if (std::fabs(area2) < 1e-9f) continue;
                pt.front = area2 < 0;           // screen y points down: outward-counter-clockwise triangles have a negative area
                if (mv.cullBack && !pt.front) continue;
                pt.sortKey = (pt.v[0].w + pt.v[1].w + pt.v[2].w) / 3.0f;
                pt.alpha = mv.alpha;
                pt.spec = mv.specular;
                dst.push_back(pt);
            }
        }
    });
    auto& target = mv.alpha < 0.999f ? transparent_ : opaque_;
    std::size_t total = 0;
    for (auto& c : out) total += c.size();
    target.reserve(target.size() + total);
    for (auto& c : out) target.insert(target.end(), c.begin(), c.end());
}

void Renderer::drawSpheres(const std::vector<Sphere>& list) {
    const bool persp = cam_.distance > 0;
    const double base = persp ? cam_.distance : 1000.0;
    const float zLimit = static_cast<float>(persp ? std::min(cam_.clipZ, base * 0.95) : cam_.clipZ);
    const auto& R = cam_.rot;
    for (const Sphere& s : list) {
        const double dx = s.x - cam_.center[0], dy = s.y - cam_.center[1], dz = s.z - cam_.center[2];
        const double vx = R[0] * dx + R[1] * dy + R[2] * dz, vy = R[3] * dx + R[4] * dy + R[5] * dz, vz = R[6] * dx + R[7] * dy + R[8] * dz;
        if (vz > zLimit) continue;
        const double w = base - vz;
        const double f = (persp ? cam_.scale * base / w : cam_.scale) * ss_;
        PreSphere p;
        p.x = static_cast<float>(W_ * 0.5 + cam_.panX * ss_ + f * vx);
        p.y = static_cast<float>(H_ * 0.5 + cam_.panY * ss_ - f * vy);
        p.w = static_cast<float>(w);
        p.oow = persp ? static_cast<float>(1.0 / w) : 1.0f;
        p.radiusPx = static_cast<float>(f * s.radius);
        p.radiusA = s.radius;
        p.r = static_cast<float>((s.color >> 16) & 255);
        p.g = static_cast<float>((s.color >> 8) & 255);
        p.b = static_cast<float>(s.color & 255);
        spheres_.push_back(p);
    }
}

void Renderer::render() {
    std::sort(transparent_.begin(), transparent_.end(), [](const PreTri& a, const PreTri& b) { return a.sortKey > b.sortKey; });   // far to near
    const int nthreads = threads_ > 0 ? threads_ : defaultThreads();
    const int bands = H_ >= 64 ? bandCount() : 1;
    parallelFor(bands, nthreads, [&](int b) {
        const int y0 = static_cast<int>(static_cast<long long>(H_) * b / bands);
        const int y1 = static_cast<int>(static_cast<long long>(H_) * (b + 1) / bands);
        rasterBand(y0, y1);
    });
    opaque_.clear();
    opaque_.shrink_to_fit();
    transparent_.clear();
    transparent_.shrink_to_fit();
    spheres_.clear();
}

void Renderer::rasterBand(int y0, int y1) const {
    const Lighting& L = light_;
    static const Normalized3 L1 = unit(-0.40f, 0.50f, 0.77f), L2 = unit(0.65f, -0.25f, 0.30f);
    static const Normalized3 H = unit(L1.x, L1.y, L1.z + 1.0f);
    const float bgR = static_cast<float>((L.background >> 16) & 255), bgG = static_cast<float>((L.background >> 8) & 255), bgB = static_cast<float>(L.background & 255);
    const bool fog = L.fogStrength > 0.0f && L.fogFar > L.fogNear;
    const float fogNear = static_cast<float>(L.fogNear), fogInv = fog ? 1.0f / static_cast<float>(L.fogFar - L.fogNear) : 0.0f;

    auto shade = [&](float nx, float ny, float nz, bool front, float r, float g, float b, float w, float specK, float& oR, float& oG, float& oB) {
        const float il = 1.0f / std::sqrt(nx * nx + ny * ny + nz * nz + 1e-20f);
        nx *= il; ny *= il; nz *= il;
        if (!front) { nx = -nx; ny = -ny; nz = -nz; }
        const float d1 = std::max(0.0f, nx * L1.x + ny * L1.y + nz * L1.z);
        const float d2 = std::max(0.0f, nx * L2.x + ny * L2.y + nz * L2.z);
        const float lit = L.ambient + L.diffuse * d1 + L.fill * d2;
        float sp = std::max(0.0f, nx * H.x + ny * H.y + nz * H.z);
        const float s2 = sp * sp, s4 = s2 * s2, s8 = s4 * s4, s16 = s8 * s8, s32 = s16 * s16;
        sp = s32 * s16 * specK * 255.0f;
        oR = r * lit + sp; oG = g * lit + sp; oB = b * lit + sp;
        if (fog) {
            const float t = std::clamp((w - fogNear) * fogInv, 0.0f, 1.0f) * L.fogStrength;
            oR += (bgR - oR) * t; oG += (bgG - oG) * t; oB += (bgB - oB) * t;
        }
    };

    auto drawTri = [&](const PreTri& t, bool blend) {
        const Vtx& a = t.v[0];
        const Vtx& b = t.v[1];
        const Vtx& c = t.v[2];
        const float minx = std::min({a.x, b.x, c.x}), maxx = std::max({a.x, b.x, c.x});
        const float miny = std::min({a.y, b.y, c.y}), maxy = std::max({a.y, b.y, c.y});
        const int ix0 = std::max(0, static_cast<int>(std::floor(minx - 0.5f)));
        const int ix1 = std::min(W_ - 1, static_cast<int>(std::ceil(maxx - 0.5f)));
        const int iy0 = std::max(y0, static_cast<int>(std::floor(miny - 0.5f)));
        const int iy1 = std::min(y1 - 1, static_cast<int>(std::ceil(maxy - 0.5f)));
        if (ix0 > ix1 || iy0 > iy1) return;
        const float area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        const float invA = 1.0f / area2;
        // lambda_i = E_i(p) / area2, E linear in p
        const float e0dx = -(c.y - b.y) * invA, e0dy = (c.x - b.x) * invA;
        const float e1dx = -(a.y - c.y) * invA, e1dy = (a.x - c.x) * invA;
        const float e2dx = -(b.y - a.y) * invA, e2dy = (b.x - a.x) * invA;
        const float px0 = ix0 + 0.5f, py0 = iy0 + 0.5f;
        float l0r = ((c.x - b.x) * (py0 - b.y) - (c.y - b.y) * (px0 - b.x)) * invA;
        float l1r = ((a.x - c.x) * (py0 - c.y) - (a.y - c.y) * (px0 - c.x)) * invA;
        float l2r = ((b.x - a.x) * (py0 - a.y) - (b.y - a.y) * (px0 - a.x)) * invA;
        for (int y = iy0; y <= iy1; y++, l0r += e0dy, l1r += e1dy, l2r += e2dy) {
            float l0 = l0r, l1 = l1r, l2 = l2r;
            std::size_t idx = static_cast<std::size_t>(y) * W_ + ix0;
            for (int x = ix0; x <= ix1; x++, idx++, l0 += e0dx, l1 += e1dx, l2 += e2dx) {
                if (l0 < 0 || l1 < 0 || l2 < 0) continue;
                const float a0 = l0 * a.oow, a1 = l1 * b.oow, a2 = l2 * c.oow;
                const float S = a0 + a1 + a2;
                const float w = (a0 * a.w + a1 * b.w + a2 * c.w) / S;
                if (w >= depth_[idx]) continue;
                const float inv = 1.0f / S;
                float r, g, bl;
                shade((a0 * a.nx + a1 * b.nx + a2 * c.nx) * inv, (a0 * a.ny + a1 * b.ny + a2 * c.ny) * inv, (a0 * a.nz + a1 * b.nz + a2 * c.nz) * inv, t.front,
                      (a0 * a.r + a1 * b.r + a2 * c.r) * inv, (a0 * a.g + a1 * b.g + a2 * c.g) * inv, (a0 * a.b + a1 * b.b + a2 * c.b) * inv, w, t.spec, r, g, bl);
                if (blend) {
                    const std::uint32_t old = color_[idx];
                    const float al = t.alpha;
                    r = r * al + static_cast<float>((old >> 16) & 255) * (1 - al);
                    g = g * al + static_cast<float>((old >> 8) & 255) * (1 - al);
                    bl = bl * al + static_cast<float>(old & 255) * (1 - al);
                } else {
                    depth_[idx] = w;
                }
                color_[idx] = pack(r, g, bl);
            }
        }
    };

    for (const PreTri& t : opaque_) {
        const float miny = std::min({t.v[0].y, t.v[1].y, t.v[2].y}), maxy = std::max({t.v[0].y, t.v[1].y, t.v[2].y});
        if (maxy < y0 - 1 || miny > y1 + 1) continue;
        drawTri(t, false);
    }

    for (const PreSphere& s : spheres_) {
        const int ix0 = std::max(0, static_cast<int>(std::floor(s.x - s.radiusPx)));
        const int ix1 = std::min(W_ - 1, static_cast<int>(std::ceil(s.x + s.radiusPx)));
        const int iy0 = std::max(y0, static_cast<int>(std::floor(s.y - s.radiusPx)));
        const int iy1 = std::min(y1 - 1, static_cast<int>(std::ceil(s.y + s.radiusPx)));
        const float invR = 1.0f / std::max(s.radiusPx, 1e-3f);
        for (int y = iy0; y <= iy1; y++) {
            const float dy = (y + 0.5f - s.y) * invR;
            for (int x = ix0; x <= ix1; x++) {
                const float dx = (x + 0.5f - s.x) * invR;
                const float r2 = dx * dx + dy * dy;
                if (r2 > 1.0f) continue;
                const float dz = std::sqrt(1.0f - r2);
                const float w = s.w - dz * s.radiusA;
                const std::size_t idx = static_cast<std::size_t>(y) * W_ + x;
                if (w >= depth_[idx]) continue;
                float r, g, b;
                shade(dx, -dy, dz, true, s.r, s.g, s.b, w, 0.35f, r, g, b);
                depth_[idx] = w;
                color_[idx] = pack(r, g, b);
            }
        }
    }

    for (const PreTri& t : transparent_) {
        const float miny = std::min({t.v[0].y, t.v[1].y, t.v[2].y}), maxy = std::max({t.v[0].y, t.v[1].y, t.v[2].y});
        if (maxy < y0 - 1 || miny > y1 + 1) continue;
        drawTri(t, true);
    }
}

std::vector<std::uint32_t> Renderer::resolve() const {
    std::vector<std::uint32_t> out(static_cast<std::size_t>(width_) * height_);
    const int n = ss_ * ss_;
    for (int y = 0; y < height_; y++)
        for (int x = 0; x < width_; x++) {
            unsigned r = 0, g = 0, b = 0;
            for (int j = 0; j < ss_; j++)
                for (int i = 0; i < ss_; i++) {
                    const std::uint32_t c = color_[static_cast<std::size_t>(y * ss_ + j) * W_ + x * ss_ + i];
                    r += (c >> 16) & 255; g += (c >> 8) & 255; b += c & 255;
                }
            out[static_cast<std::size_t>(y) * width_ + x] = 0xff000000u | ((r / n) << 16) | ((g / n) << 8) | (b / n);
        }
    return out;
}

float Renderer::depthAt(double x, double y) const {
    const int ix = static_cast<int>(x * ss_ + ss_ * 0.5), iy = static_cast<int>(y * ss_ + ss_ * 0.5);
    if (ix < 0 || iy < 0 || ix >= W_ || iy >= H_) return kFar;
    return depth_[static_cast<std::size_t>(iy) * W_ + ix];
}

}  // namespace raster
}  // namespace gui
