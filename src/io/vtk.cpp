#include "io/vtk.hpp"

#include <iomanip>

namespace ns::io {

void write_vtk_structured_grid(std::ostream& os, const StructuredMesh& mesh,
                               std::string_view title,
                               std::span<const PointScalars> point_scalars) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const long long n_nodes = static_cast<long long>(nx + 1) * (ny + 1);

    os << "# vtk DataFile Version 3.0\n";
    os << title << "\n";
    os << "ASCII\n";
    os << "DATASET STRUCTURED_GRID\n";
    os << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " 1\n";
    os << "POINTS " << n_nodes << " double\n";

    const Field<double>& x = mesh.node_x();
    const Field<double>& y = mesh.node_y();
    os << std::setprecision(16);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            os << x(i, j) << '\t' << y(i, j) << '\t' << 1.0 << '\n';
        }
    }

    if (!point_scalars.empty()) os << "POINT_DATA " << n_nodes << "\n";
    for (const auto& ps : point_scalars) {
        os << "SCALARS " << ps.name << " double\n";
        os << "LOOKUP_TABLE default\n";
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) os << (*ps.values)(i, j) << '\n';
    }
}

}  // namespace ns::io
