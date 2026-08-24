#include "mesh/mesh.hpp"

#include <charconv>
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

std::expected<StructuredMesh, Error> StructuredMesh::load_legacy(
    const std::filesystem::path& file, int nx, int ny, int ng, double scaling) {
    if (nx < 1 || ny < 1 || ng < 0)
        return std::unexpected(Error{std::format("invalid mesh dims nx={} ny={} ng={}", nx,
                                                 ny, ng)});

    TokenReader in(file);
    if (!in.ok())
        return std::unexpected(Error{std::format("grid file not readable: '{}'",
                                                 file.string())});

    const std::size_t n_nodes = static_cast<std::size_t>(nx + 1) * (ny + 1);
    StructuredMesh mesh(nx, ny, ng);

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            // skip the legacy sequential point index
            if (!in.next_double().has_value())
                return std::unexpected(Error{std::format(
                    "grid file '{}': unexpected end of data at node {} of {}", file.string(),
                    j * (nx + 1) + i, n_nodes)});
            auto x = in.next_double();
            auto y = in.next_double();
            if (!x || !y)
                return std::unexpected(Error{std::format(
                    "grid file '{}': malformed coordinate pair at node {}", file.string(),
                    j * (nx + 1) + i)});
            mesh.node_x()(i, j) = *x * scaling;
            mesh.node_y()(i, j) = *y * scaling;
        }
    }

    if (in.remaining() != 0)
        return std::unexpected(Error{std::format(
            "grid file '{}': {} trailing values after the last node", file.string(),
            in.remaining())});

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
