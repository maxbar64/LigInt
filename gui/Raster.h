// Software 3D rasterizer used by the molecular viewer: z-buffer, perspective-correct interpolation, per-pixel lighting,
// depth cueing, analytic spheres, transparent meshes, clipping and supersampling. Plain C++ (no Qt): the viewer converts the
// result to a QImage. It runs on several threads (bands of rows) and has no dependency on a GPU or on OpenGL.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "ligint/scene3d.hpp"

namespace gui {
namespace raster {

// World -> screen. The view space has x to the right, y up, z towards the viewer; the camera looks down -z.
struct Camera {
    std::array<double, 9> rot{1, 0, 0, 0, 1, 0, 0, 0, 1};   // rotation applied to (p - center), row-major
    ligint::V3 center{0, 0, 0};
    double scale = 10.0;        // pixels per Angstrom on the plane of the centre
    double panX = 0, panY = 0;  // pixels
    double distance = 0;        // camera distance from the centre (Angstrom); 0 = orthographic projection
    double clipZ = 1e9;         // view-space z above which the geometry is cut away (the front clipping plane)
    int width = 800, height = 600;

    struct Point { double x, y, w; bool ok; };   // pixel position; depth (smaller = nearer); ok = in front of the camera
    Point project(const ligint::V3& p) const;
    ligint::V3 toView(const ligint::V3& p) const;
};

// Colors are packed 0x00RRGGBB.
struct MeshView {
    const ligint::Mesh* mesh = nullptr;
    const std::vector<std::uint32_t>* colors = nullptr;   // one per vertex
    float alpha = 1.0f;                                   // < 1: transparent (drawn after the opaque geometry, sorted)
    bool cullBack = false;                                // transparent surfaces: draw only the faces that look at the viewer
    float specular = 0.25f;
};

struct Sphere { float x, y, z, radius; std::uint32_t color; };

struct Lighting {
    float ambient = 0.30f, diffuse = 0.62f, fill = 0.14f;
    std::uint32_t background = 0xffffff;
    float fogStrength = 0.0f;                             // 0..1: fraction of the background mixed in at the far end
    double fogNear = 0, fogFar = 1;                       // depth range of the fog (same units as Point::w)
};

class Renderer {
public:
    // supersample: 1 = none, 2 = 2x2 samples per pixel, ...
    Renderer(int width, int height, int supersample = 1);

    void setCamera(const Camera& camera);
    void setLighting(const Lighting& lighting);          // also clears the image to the background: call it before drawing
    void setThreads(int n) { threads_ = n; }              // 0 = automatic

    // Draw calls only queue the work (the transformation of the vertices runs at once, on several threads).
    void drawMesh(const MeshView& mesh);
    void drawSpheres(const std::vector<Sphere>& spheres);

    // Rasterizes everything that was queued: opaque meshes and spheres first, then the transparent meshes.
    void render();

    // 0xFFRRGGBB, width * height pixels, after averaging the samples.
    std::vector<std::uint32_t> resolve() const;
    // Depth (Camera::Point::w) at an output pixel; a very large value where nothing was drawn.
    float depthAt(double x, double y) const;

    int width() const { return width_; }
    int height() const { return height_; }

private:
    struct Vtx { float x, y, oow, w, nx, ny, nz, r, g, b; };
    struct PreTri { Vtx v[3]; bool front; float sortKey; float alpha; float spec; };
    struct PreSphere { float x, y, oow, w, radiusPx, radiusA, r, g, b; };

    void rasterBand(int y0, int y1) const;
    int bandCount() const;

    int width_, height_, ss_, W_, H_, threads_ = 0;
    Camera cam_;
    Lighting light_;
    std::vector<PreTri> opaque_, transparent_;
    std::vector<PreSphere> spheres_;
    mutable std::vector<float> depth_;
    mutable std::vector<std::uint32_t> color_;   // 0x00RRGGBB per sample
};

// Fills `out` (width*height ARGB) with the image; convenience used by the tests.
inline std::uint32_t pack(float r, float g, float b) {
    auto c = [](float v) { return static_cast<std::uint32_t>(v < 0 ? 0 : (v > 255 ? 255 : v + 0.5f)); };
    return (c(r) << 16) | (c(g) << 8) | c(b);
}

}  // namespace raster
}  // namespace gui
