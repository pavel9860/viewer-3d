// Geometry shared by the DXF entities: object coordinate systems, curve sampling, NURBS,
// polygon triangulation.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace cadcore::formats::autocad {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    bool operator==(const Vec3&) const = default;
};
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double length(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalized(const Vec3& a) { double l = length(a); return l > 0 ? a * (1.0 / l) : Vec3{0, 0, 1}; }

// Object coordinate system of an extrusion direction (DXF "arbitrary axis algorithm").
struct Ocs {
    Vec3 ax{1, 0, 0}, ay{0, 1, 0}, az{0, 0, 1};
    explicit Ocs(Vec3 extrusion = {0, 0, 1});
    bool identity() const { return az == Vec3{0, 0, 1} && ax == Vec3{1, 0, 0}; }
    Vec3 to_wcs(const Vec3& p) const { return ax * p.x + ay * p.y + az * p.z; }
};

using Polyline3 = std::vector<Vec3>;

// Segments for an arc of `sweep` radians and radius r, with chord deviation at most
// kChordTolerance * r (at least 1 segment; a full circle gets 72).
int arc_segments(double sweep);
constexpr double kChordTolerance = 1e-3;

// Points of the arc centre c (2D, z = elevation), radius r, from angle a0 sweeping `sweep`
// radians (positive = counter-clockwise); includes both ends.
void append_arc(Polyline3& out, Vec3 c, double r, double a0, double sweep, bool skip_first);

// Points of the arc a polyline bulge b describes from p0 to p1 (b = tan(angle / 4), positive =
// counter-clockwise), without p0, with p1. A zero bulge appends p1 only.
void append_bulge(Polyline3& out, Vec3 p0, Vec3 p1, double bulge);

// Points of the ellipse arc centre c, semi-axes u, v (vectors), parameters t0..t1 (t1 > t0).
void append_ellipse(Polyline3& out, Vec3 c, Vec3 u, Vec3 v, double t0, double t1, bool skip_first);

// Points of a NURBS curve (degree, knots, control points, weights empty = 1), sampled per
// non-empty knot span. Returns false when the data is inconsistent.
bool sample_nurbs(Polyline3& out, int degree, const std::vector<double>& knots, const std::vector<Vec3>& ctrl,
                  const std::vector<double>& weights);

// Control points and knots of the cubic (or lower) B-spline interpolating `fit` points with
// chord-length parameters; end tangents optional (zero vector = not given).
bool interpolate_spline(const std::vector<Vec3>& fit, int degree, Vec3 t0, Vec3 t1, std::vector<double>& knots,
                        std::vector<Vec3>& ctrl);

// Fills the region bounded by closed loops in the z = 0 plane under the even-odd rule (a point is
// inside when a ray from it crosses the loops an odd number of times). Loops may nest, touch and
// cross. Appends triangles (3 points each, counter-clockwise, z = 0): the region is cut into
// horizontal slabs at every vertex and edge crossing; within a slab the edges do not cross, so
// the region there is a set of trapezoids.
void fill_even_odd(const std::vector<Polyline3>& loops, std::vector<Vec3>& tris);

// Nesting depth of each loop: how many other loops contain its first vertex.
std::vector<int> nesting_depths(const std::vector<Polyline3>& loops);

}  // namespace cadcore::formats::autocad
