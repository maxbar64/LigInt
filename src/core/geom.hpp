// Small vector helpers shared by the geometry code (scene3d). Internal, not installed.
#pragma once

#include <array>
#include <cmath>

#include "ligint/ligint.hpp"

namespace ligint {
namespace geom {

using V3 = ligint::V3;

inline V3 add(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline V3 mul(const V3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline V3 cross(const V3& a, const V3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double len(const V3& a) { return std::sqrt(dot(a, a)); }
inline double dist(const V3& a, const V3& b) { return len(sub(a, b)); }
inline V3 normalized(const V3& a, const V3& fallback = {0, 0, 1}) {
    const double l = len(a);
    return l > 1e-12 ? mul(a, 1.0 / l) : fallback;
}
inline V3 lerp(const V3& a, const V3& b, double t) { return add(mul(a, 1.0 - t), mul(b, t)); }

// Any unit vector perpendicular to v.
inline V3 perpendicular(const V3& v) {
    const V3 axis = std::fabs(v[0]) < 0.6 ? V3{1, 0, 0} : V3{0, 1, 0};
    return normalized(cross(v, axis));
}

}  // namespace geom
}  // namespace ligint
