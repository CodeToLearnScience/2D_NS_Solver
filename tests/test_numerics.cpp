// Parity tests for the inviscid scheme pipelines: complete LLF chains must
// reproduce legacy outputs bit-for-bit. Goldens in goldens_numerics.inc were
// generated from the legacy objects (see REFACTORING_PLAN.md, parity
// methodology).

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "numerics/face_data.hpp"
#include "numerics/llf.hpp"

namespace {

using namespace ns;
using namespace ns::numerics;
using namespace ns::physics;
#include "goldens_numerics.inc"

constexpr int kNx = 60, kNy = 20, kNg = 2;
constexpr const char* kGridPath = NS_SOURCE_DIR "/input/RampGrid6020.dat";

struct SchemeRig {
    StructuredMesh mesh{kNx, kNy, kNg};
    MeshMetrics metrics{compute_metrics(mesh)};
    MultiField<double> dv{8, kNx, kNy, kNg};
    IdealGas eos{};
    FaceFluxData ifaces{kNConv, kNx + 1, kNy + 1, kNg};
    FaceFluxData jfaces{kNConv, kNx + 1, kNy + 1, kNg};

    SchemeRig() {
        for (int i = -kNg; i < kNx + kNg; ++i)
            for (int j = -kNg; j < kNy + kNg; ++j) {
                const double a = 0.3 * (i + 2), b = 0.47 * (j + 2);
                const double w1 = std::sin(a) * std::cos(b);
                const double w2 = std::cos(0.21 * a) * std::sin(0.53 * b);
                dv(0, i, j) = 1.0 + 0.2 * w1;
                dv(1, i, j) = 120.0 * (0.5 + 0.3 * w2);
                dv(2, i, j) = -35.0 + 12.0 * w1;
                dv(4, i, j) = 260.0 + 80.0 * w1;
                dv(3, i, j) = 101325.0 * (0.8 + 0.05 * w2);
                dv(5, i, j) = 310.0 + 15.0 * w2;
                dv(6, i, j) = 1.8e-5 * (1.0 + 0.4 * w2);
                dv(7, i, j) = 0.9 * (1.0 + 0.3 * w1);
            }
        define_convective_faces(mesh, metrics, dv, eos, ifaces, jfaces);
    }

    void reload_grid() {
        auto loaded = StructuredMesh::load_legacy(kGridPath, kNx, kNy, kNg);
        if (!loaded) throw std::runtime_error(loaded.error().message);
        mesh = std::move(*loaded);
        metrics = compute_metrics(mesh);
        define_convective_faces(mesh, metrics, dv, eos, ifaces, jfaces);
    }
};

TEST(NumericsParity, LlfFirstOrderMatchesLegacyBitwise) {
    SchemeRig rig;
    rig.reload_grid();

    MultiField<double> diss(kNConv, kNx, kNy, kNg);
    diss.fill(0.0);
    llf_dissipation_first_order(rig.mesh, rig.metrics, rig.dv, rig.ifaces.con_var_diff,
                                rig.jfaces.con_var_diff, rig.eos, diss);

    for (const auto& g : goldens::diss_llf1)
        for (int k = 0; k < 4; ++k)
            EXPECT_DOUBLE_EQ(diss(k, g.i, g.j), g.v[k])
                << "diss_llf1 plane " << k << " at (" << g.i << "," << g.j << ")";

    MultiField<double> rhs(kNConv, kNx, kNy, kNg);
    assemble_rhs_from_average_fluxes(rig.mesh, rig.metrics, rig.ifaces.avg_flux,
                                     rig.jfaces.avg_flux, diss, rhs);
    for (const auto& g : goldens::rhs_llf1)
        for (int k = 0; k < 4; ++k)
            EXPECT_DOUBLE_EQ(rhs(k, g.i, g.j), g.v[k])
                << "rhs_llf1 plane " << k << " at (" << g.i << "," << g.j << ")";
}

TEST(NumericsParity, LlfSecondOrderMatchesLegacyBitwise) {
    SchemeRig rig;
    rig.reload_grid();

    MultiField<double> dui(kNConv, kNx, kNy, kNg), duj(kNConv, kNx, kNy, kNg);
    primitive_differences(rig.dv, dui, duj);

    MultiField<double> diss(kNConv, kNx, kNy, kNg);
    diss.fill(0.0);
    llf_dissipation_second_order(rig.mesh, rig.metrics, rig.dv, dui, duj, rig.eos, diss);

    for (const auto& g : goldens::diss_llf2)
        for (int k = 0; k < 4; ++k)
            EXPECT_DOUBLE_EQ(diss(k, g.i, g.j), g.v[k])
                << "diss_llf2 plane " << k << " at (" << g.i << "," << g.j << ")";

    MultiField<double> rhs(kNConv, kNx, kNy, kNg);
    assemble_rhs_from_average_fluxes(rig.mesh, rig.metrics, rig.ifaces.avg_flux,
                                     rig.jfaces.avg_flux, diss, rhs);
    for (const auto& g : goldens::rhs_llf2)
        for (int k = 0; k < 4; ++k)
            EXPECT_DOUBLE_EQ(rhs(k, g.i, g.j), g.v[k])
                << "rhs_llf2 plane " << k << " at (" << g.i << "," << g.j << ")";
}

// ---------------------------------------------------------------------------
// Invariant tests
// ---------------------------------------------------------------------------
TEST(NumericsInvariants, ConstantFieldGivesZeroRhsBothOrders) {
    StructuredMesh mesh{8, 6, 3};
    for (int i = -3; i <= 9; ++i)
        for (int j = -3; j <= 7; ++j) {
            mesh.node_x()(i, j) = i * 0.5;
            mesh.node_y()(i, j) = j * 0.5 + 0.01 * i;  // mild skew
        }
    MeshMetrics m = compute_metrics(mesh);

    MultiField<double> dv(8, 8, 6, 3);
    for (int i = -3; i < 11; ++i)
        for (int j = -3; j < 9; ++j) {
            dv(physics::dvi::RHO, i, j) = 1.1;
            dv(physics::dvi::U, i, j) = 250.0;
            dv(physics::dvi::V, i, j) = 40.0;
            dv(physics::dvi::P, i, j) = 80000.0;
        }

    IdealGas eos{};
    FaceFluxData fi(kNConv, 9, 7, 3), fj(kNConv, 9, 7, 3);
    define_convective_faces(mesh, m, dv, eos, fi, fj);

    // Flux consistency: F(U)*n equal on both sides of every face.
    for (int f = 0; f <= 8; ++f)
        for (int j = 0; j < 6; ++j)
            for (int k = 0; k < 4; ++k) EXPECT_DOUBLE_EQ(fi.flux_diff(k, f, j), 0.0);

    MultiField<double> diss(kNConv, 8, 6, 3);
    diss.fill(0.0);
    llf_dissipation_first_order(mesh, m, dv, fi.con_var_diff, fj.con_var_diff, eos, diss);
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 6; ++j) EXPECT_DOUBLE_EQ(diss(k, i, j), 0.0);

    MultiField<double> rhs(kNConv, 8, 6, 3);
    assemble_rhs_from_average_fluxes(mesh, m, fi.avg_flux, fj.avg_flux, diss, rhs);

    // Roundoff of the telescoping sum scales with the largest scattered term
    // (energy fluxes are O(1e6) here), so bound residuals relatively.
    double scale = 0.0;
    for (int k = 0; k < 4; ++k)
        for (int f = 0; f <= 8; ++f)
            for (int j = 0; j < 6; ++j)
                scale = std::max({scale, std::fabs(fi.avg_flux(k, f, j)),
                                  std::fabs(fj.avg_flux(k, f, j))});
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 6; ++j)
                // The vcont average-state flux is not exactly flux-consistent
                // under skewing; free-stream preservation holds to roundoff.
                EXPECT_NEAR(rhs(k, i, j), 0.0, 1e-13 * scale);
}

TEST(NumericsInvariants, MinmodFluxLimiterBoundsAndValues) {
    EXPECT_DOUBLE_EQ(minmod_flux_limiter(1.0, 1.0), 1.0);
    EXPECT_DOUBLE_EQ(minmod_flux_limiter(-1.0, 1.0), 0.0);   // sign change -> 0
    EXPECT_DOUBLE_EQ(minmod_flux_limiter(2.0, 1.0), 1.0);    // r > 1 clamps
    EXPECT_NEAR(minmod_flux_limiter(0.5, 1.0), 0.5 / (1.0 + 1e-16), 1e-15);
}

}  // namespace
