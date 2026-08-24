#include "mesh/mesh.hpp"

#include <algorithm>
#include <charconv>
#include <utility>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>

namespace ns {
namespace {

class TokenReader {
public:
    explicit TokenReader(const std::filesystem::path& p) {
        std::ifstream in(p);
        if (!in) return;
        std::ostringstream raw;
        raw << in.rdbuf();
        std::istringstream lines(raw.str());
        std::string line, tok;
        while (std::getline(lines, line)) {
            if (auto pos = line.find('#'); pos != std::string::npos) line.resize(pos);
            std::istringstream ls(line);
            while (ls >> tok) tokens_.push_back(tok);
        }
        ok_ = true;
    }

    [[nodiscard]] bool ok() const { return ok_; }
    [[nodiscard]] std::size_t remaining() const { return tokens_.size() - pos_; }

    std::optional<double> next_double() {
        if (pos_ >= tokens_.size()) return std::nullopt;
        const std::string& s = tokens_[pos_++];
        double v{};
        auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
        if (ec != std::errc{}) return std::nullopt;
        return v;
    }

private:
    bool ok_ = false;
    std::vector<std::string> tokens_;
    std::size_t pos_ = 0;
};

}  // namespace

StructuredMesh::StructuredMesh(int nx, int ny, int ng)
    : nx_(nx), ny_(ny), ng_(ng), x_(nx + 1, ny + 1, ng), y_(nx + 1, ny + 1, ng) {}

namespace {
// Reads every node pair of the file into flat arrays (x,y interleaved per node).
std::expected<std::pair<std::vector<double>, std::vector<double>>, Error> read_all_nodes(
    const std::filesystem::path& file, int nx, int ny) {
    TokenReader in(file);
    if (!in.ok())
        return std::unexpected(
            Error{std::format("grid file not readable: '{}'", file.string())});
    std::vector<double> xs, ys;
    const std::size_t n =
        static_cast<std::size_t>(nx + 1) * static_cast<std::size_t>(ny + 1);
    xs.reserve(n);
    ys.reserve(n);
    for (std::size_t p = 0; p < n; ++p) {
        if (!in.next_double().has_value())
            return std::unexpected(Error{std::format(
                "grid file '{}': unexpected end of data at node {} of {}", file.string(),
                p, n)});
        auto x = in.next_double();
        auto y = in.next_double();
        if (!x || !y)
            return std::unexpected(Error{std::format(
                "grid file '{}': malformed coordinate pair at node {}", file.string(), p)});
        xs.push_back(*x);
        ys.push_back(*y);
    }
    if (in.remaining() != 0)
        return std::unexpected(Error{std::format(
            "grid file '{}': {} trailing values after the last node", file.string(),
            in.remaining())});
    return std::make_pair(std::move(xs), std::move(ys));
}
}  // namespace

std::expected<StructuredMesh, Error> StructuredMesh::load_legacy(
    const std::filesystem::path& file, int nx, int ny, int ng, double scaling) {
    if (nx < 1 || ny < 1 || ng < 0)
        return std::unexpected(Error{std::format("invalid mesh dims nx={} ny={} ng={}", nx,
                                                 ny, ng)});
    auto nodes = read_all_nodes(file, nx, ny);
    if (!nodes) return std::unexpected(nodes.error());
    const auto& xs = nodes->first;
    const auto& ys = nodes->second;

    StructuredMesh mesh(nx, ny, ng);
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            const std::size_t p = static_cast<std::size_t>(j) * (nx + 1) + i;
            mesh.node_x()(i, j) = xs[p] * scaling;
            mesh.node_y()(i, j) = ys[p] * scaling;
        }
    return mesh;
}

std::expected<StructuredMesh, Error> StructuredMesh::load_legacy_slab(
    const std::filesystem::path& file, int nx, int ny_global, int j0, int j1, int ng,
    double scaling) {
    if (nx < 1 || ny_global < 1 || ng < 0 || j0 < 0 || j1 <= j0 || j1 > ny_global)
        return std::unexpected(
            Error{std::format("invalid slab request nx={} ny={} [{},{}) ng={}", nx,
                              ny_global, j0, j1, ng)});
    auto nodes = read_all_nodes(file, nx, ny_global);
    if (!nodes) return std::unexpected(nodes.error());
    const auto& xs = nodes->first;
    const auto& ys = nodes->second;

    const int ny_local = j1 - j0;
    StructuredMesh mesh(nx, ny_local, ng);
    for (int jl = -ng; jl <= ny_local + ng; ++jl) {
        const int jg = std::clamp(j0 + jl, 0, ny_global);
        for (int il = -ng; il <= nx + ng; ++il) {
            const int ig = std::clamp(il, 0, nx);
            const std::size_t p = static_cast<std::size_t>(jg) * (nx + 1) + ig;
            mesh.node_x()(il, jl) = xs[p] * scaling;
            mesh.node_y()(il, jl) = ys[p] * scaling;
        }
    }
    return mesh;
}

MeshMetrics compute_metrics(const StructuredMesh& mesh) {
    const int nx = mesh.nx(), ny = mesh.ny(), ng = mesh.ng();
    MeshMetrics m{
        .si = Field<Vec2>(nx, ny, ng),
        .sj = Field<Vec2>(nx, ny, ng),
        .area = Field<double>(nx, ny, ng),
    };
    const Field<double>& x = mesh.node_x();
    const Field<double>& y = mesh.node_y();

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const double dxj = x(i, j + 1) - x(i, j);   // tangent along +j at the left face
            const double dyj = y(i, j + 1) - y(i, j);
            m.si(i, j) = {-dyj, dxj};

            const double dxi = x(i + 1, j) - x(i, j);   // tangent along +i at the bottom face
            const double dyi = y(i + 1, j) - y(i, j);
            m.sj(i, j) = {dyi, -dxi};

            m.area(i, j) = 0.5 * ((x(i, j) - x(i + 1, j + 1)) * (y(i + 1, j) - y(i, j + 1)) +
                                  (x(i, j + 1) - x(i + 1, j)) * (y(i, j) - y(i + 1, j + 1)));
        }
    }

    m.si.fill_ghosts_copy();
    m.sj.fill_ghosts_copy();
    m.area.fill_ghosts_copy();
    return m;
}

}  // namespace ns
