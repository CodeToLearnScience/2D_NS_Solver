// Unit tests for StructuredMesh: legacy grid reader and geometric metrics,
// verified against hand-computed geometry.

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "mesh/mesh.hpp"

namespace {

namespace fs = std::filesystem;

int g_counter = 0;

fs::path write_grid_file(const std::string& content) {
    const auto path = fs::temp_directory_path() /
                      std::format("ns_mesh_{}_{}.dat", ::getpid(), g_counter++);
    std::ofstream out(path);
    out << content;
    return path;
}

TEST(MeshLoad, ReadsLegacyIxFastestFormat) {
    // 3x2 cells -> 4x3 nodes; i fastest per row. Comments must be skipped.
    const auto path = write_grid_file(
        "# comment line\n"
        "1 0.0 0.0\n2 1.0 0.0\n3 2.0 0.0\n4 3.0 0.0\n"
        "\n"
        "5 0.0 1.0\n6 1.5 1.0\n7 2.0 2.5\n8 3.0 1.0\n"   // note skewed middle node
        "9 0.0 2.0\n10 1.0 2.0\n11 2.0 2.0\n12 3.0 2.0\n");

    auto mesh = ns::StructuredMesh::load_legacy(path, 3, 2, 2);
    ASSERT_TRUE(mesh.has_value()) << mesh.error().message;

    EXPECT_DOUBLE_EQ(mesh->node_x()(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(mesh->node_x()(2, 0), 2.0);   // x(2,0) = third value of row j=0
    EXPECT_DOUBLE_EQ(mesh->node_y()(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(mesh->node_x()(1, 1), 1.5);
    EXPECT_DOUBLE_EQ(mesh->node_y()(1, 1 + 1), 2.0);
}

TEST(MeshLoad, AppliesScaling) {
    const auto path = write_grid_file(
        "1 0 0\n2 2 0\n3 4 0\n4 6 0\n"
        "5 0 2\n6 2 2\n7 4 2\n8 6 2\n");
    auto mesh = ns::StructuredMesh::load_legacy(path, 3, 1, 1, 0.5);
    ASSERT_TRUE(mesh.has_value());
    EXPECT_DOUBLE_EQ(mesh->node_x()(3, 0), 3.0);
    EXPECT_DOUBLE_EQ(mesh->node_y()(0, 1), 1.0);
}

TEST(MeshLoad, RejectsTruncatedAndOverlongFiles) {
    auto short_path = write_grid_file("1 0 0\n2 1 0\n");
    auto bad = ns::StructuredMesh::load_legacy(short_path, 3, 1, 1);
    ASSERT_FALSE(bad.has_value());
    EXPECT_NE(bad.error().message.find("unexpected end"), std::string::npos);

    auto long_path = write_grid_file(
        "1 0 0\n2 1 0\n3 0 1\n4 1 1\n5 9 9\n");  // extra point
    bad = ns::StructuredMesh::load_legacy(long_path, 1, 1, 1);
    ASSERT_FALSE(bad.has_value());
    EXPECT_NE(bad.error().message.find("trailing values"), std::string::npos);
}

TEST(MeshMetrics, UniformCartesianGrid) {
    // 4x3 unit-spaced cells: left-face vector of every cell is (-1, 0),
    // bottom-face vector is (0, -1) (both outward), area is 1.
    auto mesh = ns::StructuredMesh(4, 3, 2);
    for (int i = 0; i <= 4; ++i)
        for (int j = 0; j <= 3; ++j) {
            mesh.node_x()(i, j) = i;
            mesh.node_y()(i, j) = j;
        }

    const auto m = ns::compute_metrics(mesh);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(m.si(i, j).x, -1.0, 1e-14);
            EXPECT_NEAR(m.si(i, j).y, 0.0, 1e-14);
            EXPECT_NEAR(m.sj(i, j).x, 0.0, 1e-14);
            EXPECT_NEAR(m.sj(i, j).y, -1.0, 1e-14);
            EXPECT_NEAR(m.area(i, j), 1.0, 1e-14);
        }
    }
    EXPECT_DOUBLE_EQ(m.area(-1, -1), m.area(0, 0));  // ghost copy
}

TEST(MeshMetrics, AffineCellAreaEqualsJacobianDet) {
    // u = a*xi + b*eta ; v = c*xi + d*eta  -> every cell has |area| = |ad-bc|,
    // face vectors are constant: si = (-b, d)*h? derived directly below.
    const double h = 0.25;
    auto mesh = ns::StructuredMesh(2, 2, 1);
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
            mesh.node_x()(i, j) = 2.0 * i * h + 0.7 * j * h;
            mesh.node_y()(i, j) = -0.3 * i * h + 1.4 * j * h;
        }
    const double det = (2.0 * 1.4 - 0.7 * (-0.3)) * h * h;  // = 3.01 h^2

    const auto m = ns::compute_metrics(mesh);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            EXPECT_NEAR(m.area(i, j), det, 1e-13);

    // Left-face tangent along +j is (0.7, 1.4)*h -> si = (-1.4, 0.7)*h.
    EXPECT_NEAR(m.si(1, 1).x, -1.4 * h, 1e-13);
    EXPECT_NEAR(m.si(1, 1).y, 0.7 * h, 1e-13);
    // Bottom-face tangent along +i is (2.0, -0.3)*h -> sj = (-0.3, -2.0)*h.
    EXPECT_NEAR(m.sj(1, 1).x, -0.3 * h, 1e-13);
    EXPECT_NEAR(m.sj(1, 1).y, -2.0 * h, 1e-13);
}

TEST(MeshMetrics, GeneralQuadMatchesShoelace) {
    // A single badly-skewed cell: metrics area must equal the shoelace area.
    auto mesh = ns::StructuredMesh(1, 1, 1);
    mesh.node_x()(0, 0) = 0.0; mesh.node_y()(0, 0) = 0.0;
    mesh.node_x()(1, 0) = 2.5; mesh.node_y()(1, 0) = 0.25;
    mesh.node_x()(0, 1) = -0.5; mesh.node_y()(0, 1) = 1.75;
    mesh.node_x()(1, 1) = 1.75; mesh.node_y()(1, 1) = 2.4;

    const auto m = ns::compute_metrics(mesh);
    // Shoelace over A(0,0) -> B(2.5,0.25) -> C(1.75,2.4) -> D(-0.5,1.75).
    const double shoelace =
        0.5 * ((0.0 * 0.25 - 2.5 * 0.0) +          // A x B
               (2.5 * 2.4 - 1.75 * 0.25) +         // B x C
               (1.75 * 1.75 - (-0.5) * 2.4) +      // C x D
               ((-0.5) * 0.0 - 0.0 * 1.75));       // D x A
    EXPECT_NEAR(m.area(0, 0), shoelace, 1e-12);
    EXPECT_GT(m.area(0, 0), 0.0);
}

TEST(MeshRealData, HalfCylinder45LoadsAndHasPositiveArea) {
    const auto grid = NS_SOURCE_DIR "/input/HalfCylinderGrid4545.dat";
    if (!fs::exists(grid)) GTEST_SKIP() << "grid fixture not present";

    auto loaded = ns::StructuredMesh::load_legacy(grid, 45, 45, 3);
    ASSERT_TRUE(loaded.has_value()) << loaded.error().message;
    ns::StructuredMesh& mesh = *loaded;

    const auto m = ns::compute_metrics(mesh);
    double total = 0.0;
    for (int i = 0; i < 45; ++i)
        for (int j = 0; j < 45; ++j) {
            ASSERT_GT(m.area(i, j), 0.0) << "non-positive area at " << i << "," << j;
            total += m.area(i, j);
        }
    EXPECT_GT(total, 0.0);

    // First file node was (~0, 1): x ~ 6.12e-17, y = 1.
    EXPECT_NEAR(mesh.node_x()(0, 0), 6.12323399573677e-17, 1e-30);
    EXPECT_DOUBLE_EQ(mesh.node_y()(0, 0), 1.0);
}

}  // namespace
