// Unit tests for the IO layer: VTK structured-grid writer and binary restart.

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

#include "io/restart.hpp"
#include "io/vtk.hpp"

namespace {

namespace fs = std::filesystem;

ns::StructuredMesh make_test_mesh() {
    ns::StructuredMesh mesh(2, 1, 2);  // 2x1 cells, 3x2 nodes
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 1; ++j) {
            mesh.node_x()(i, j) = i;
            mesh.node_y()(i, j) = j;
        }
    return mesh;
}

TEST(VtkWriter, EmitsLegacyCompatibleStructuredGrid) {
    auto mesh = make_test_mesh();
    ns::Field<double> rho(mesh.nx() + 1, mesh.ny() + 1, mesh.ng());
    rho.fill(1.25);
    const std::string title = "unit-test field";

    const std::array scalars{ns::io::PointScalars{std::string("density"), &rho}};

    std::ostringstream os;
    ns::io::write_vtk_structured_grid(os, mesh, title, scalars);
    const std::string text = os.str();

    EXPECT_NE(text.find("# vtk DataFile Version 3.0"), std::string::npos);
    EXPECT_NE(text.find(title), std::string::npos);
    EXPECT_NE(text.find("DATASET STRUCTURED_GRID"), std::string::npos);
    EXPECT_NE(text.find("DIMENSIONS 3 2 1"), std::string::npos);
    EXPECT_NE(text.find("POINTS 6 double"), std::string::npos);
    // Node ordering j-outer/i-inner: first line is node (0,0), second is (1,0).
    EXPECT_NE(text.find("\n0\t0\t1\n1\t0\t1\n"), std::string::npos);
    EXPECT_NE(text.find("POINT_DATA 6"), std::string::npos);
    EXPECT_NE(text.find("SCALARS density double"), std::string::npos);
    EXPECT_NE(text.find("LOOKUP_TABLE default"), std::string::npos);
}

TEST(RestartIo, RoundTripIsBitExact) {
    auto mesh = make_test_mesh();
    ns::Field<double> a(mesh.nx(), mesh.ny(), mesh.ng());
    ns::Field<double> b(mesh.nx(), mesh.ny(), mesh.ng());
    for (int i = -mesh.ng(); i < mesh.nx() + mesh.ng(); ++i)
        for (int j = -mesh.ng(); j < mesh.ny() + mesh.ng(); ++j) {
            a(i, j) = std::sin(3.0 * i + 1.0) * std::cos(j * 2.0);
            b(i, j) = i * 100.0 + j;
        }

    const auto path =
        fs::temp_directory_path() / std::format("ns_restart_{}.bin", ::getpid());
    const std::vector<std::pair<std::string, const ns::Field<double>*>> out_fields{
        {std::string("cv_rho"), &a}, {std::string("x-velocity"), &b}};
    const auto wrote = ns::io::write_restart(path, mesh.nx(), mesh.ny(), mesh.ng(),
                                             out_fields);
    ASSERT_TRUE(wrote.has_value()) << wrote.error().message;

    const auto loaded = ns::io::read_restart(path);
    ASSERT_TRUE(loaded.has_value()) << loaded.error().message;
    EXPECT_EQ(loaded->nx, 2);
    EXPECT_EQ(loaded->ny, 1);
    EXPECT_EQ(loaded->ng, 2);
    ASSERT_EQ(loaded->fields.size(), 2u);

    const ns::io::RestartField* ra = loaded->find("cv_rho");
    ASSERT_NE(ra, nullptr);
    ASSERT_EQ(ra->values.size(), a.size());
    EXPECT_EQ(std::memcmp(ra->values.data(), a.data(), a.size() * sizeof(double)), 0);

    const ns::io::RestartField* rb = loaded->find("x-velocity");
    ASSERT_NE(rb, nullptr);
    EXPECT_EQ(std::memcmp(rb->values.data(), b.data(), b.size() * sizeof(double)), 0);

    EXPECT_EQ(loaded->find("missing"), nullptr);
}

TEST(RestartIo, RejectsNonRestartFiles) {
    const auto path =
        fs::temp_directory_path() / std::format("ns_restart_bad_{}.bin", ::getpid());
    { std::ofstream(path) << "definitely not a restart file"; }
    const auto res = ns::io::read_restart(path);
    ASSERT_FALSE(res.has_value());
    EXPECT_NE(res.error().message.find("not a solver restart file"), std::string::npos);
}

TEST(RestartIo, RejectsTruncatedPayload) {
    auto mesh = make_test_mesh();
    ns::Field<double> a(mesh.nx(), mesh.ny(), mesh.ng());
    a.fill(1.0);
    const auto path =
        fs::temp_directory_path() / std::format("ns_restart_trunc_{}.bin", ::getpid());
    const std::vector<std::pair<std::string, const ns::Field<double>*>> one{
        {std::string("f"), &a}};
    ASSERT_TRUE(
        ns::io::write_restart(path, mesh.nx(), mesh.ny(), mesh.ng(), one).has_value());

    // Corrupt by truncating half the payload.
    const auto truncated = fs::temp_directory_path() /
                           std::format("ns_restart_trunc2_{}.bin", ::getpid());
    {
        std::ifstream in(path, std::ios::binary);
        std::string data((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
        in.close();
        std::ofstream out(truncated, std::ios::binary);
        out.write(data.data(), static_cast<std::streamsize>(data.size() / 2));
    }
    const auto res = ns::io::read_restart(truncated);
    ASSERT_FALSE(res.has_value());
}

}  // namespace
