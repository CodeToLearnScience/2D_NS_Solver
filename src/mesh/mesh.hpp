#pragma once

// Structured body-fitted mesh: node coordinates, legacy-format grid reader,
// and geometric metrics (face vectors + cell areas) matching the legacy
// formulas exactly (see src/mesh/metrics in grid_computations.cpp):
//
//   si(i,j) = { y(i,j) - y(i,j+1),  x(i,j+1) - x(i,j) }   // left i-face
//   sj(i,j) = { y(i+1,j) - y(i,j),  x(i,j) - x(i+1,j) }   // bottom j-face
//
// These are outward-pointing face vectors; |si| is the face length.

#include <expected>
#include <filesystem>
#include <string>

#include "core/error.hpp"
#include "fields/field.hpp"

namespace ns {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

class StructuredMesh {
public:
    StructuredMesh(int nx, int ny, int ng);

    /// Reads the legacy grid format: (nx+1)*(ny+1) whitespace-separated
    /// "index x y" triplets, i fastest within each j row, optional '#' comments.
    static std::expected<StructuredMesh, Error> load_legacy(
        const std::filesystem::path& file, int nx, int ny, int ng, double scaling = 1.0);

    /// MPI helper: reads the full grid file but keeps only cell rows
    /// [j0, j1) as the local interior; node rows outside
    /// [max(0,j0-ng), min(ny_global-1, j1-1+ng)+1] are clamp-filled from the
    /// nearest global row (ghost node coordinates never feed the kernels).
    static std::expected<StructuredMesh, Error> load_legacy_slab(
        const std::filesystem::path& file, int nx, int ny_global, int j0, int j1,
        int ng, double scaling = 1.0);

    [[nodiscard]] int nx() const noexcept { return nx_; }
    [[nodiscard]] int ny() const noexcept { return ny_; }
    [[nodiscard]] int ng() const noexcept { return ng_; }

    /// Node coordinates on an (nx+1) x (ny+1) interior canvas (ghosted).
    [[nodiscard]] const Field<double>& node_x() const noexcept { return x_; }
    [[nodiscard]] const Field<double>& node_y() const noexcept { return y_; }
    Field<double>& node_x() noexcept { return x_; }
    Field<double>& node_y() noexcept { return y_; }

private:
    int nx_, ny_, ng_;
    Field<double> x_, y_;
};

struct MeshMetrics {
    Field<Vec2> si;        // i-face vector at the LEFT face of each cell
    Field<Vec2> sj;        // j-face vector at the BOTTOM face of each cell
    Field<double> area;    // signed cell area (positive for right-handed ij)
};

/// Computes metrics on physical cells and fills ghost layers by copy.
[[nodiscard]] MeshMetrics compute_metrics(const StructuredMesh& mesh);

}  // namespace ns
