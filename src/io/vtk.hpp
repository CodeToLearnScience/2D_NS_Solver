#pragma once

// ASCII VTK structured-grid writer (legacy-compatible layout:
// DATASET STRUCTURED_GRID, j-outer/i-inner node ordering).

#include <ostream>
#include <span>
#include <string>
#include <string_view>

#include "fields/field.hpp"
#include "mesh/mesh.hpp"

namespace ns::io {

struct PointScalars {
    std::string name;
    const Field<double>* values = nullptr;  // node-canvas data
};

void write_vtk_structured_grid(std::ostream& os, const StructuredMesh& mesh,
                               std::string_view title,
                               std::span<const PointScalars> point_scalars);

}  // namespace ns::io
