// Round-trip tests: legacy Input*.cfg + GridTop*.dat pairs from the repo must
// import into the typed Config model with values matching the original files.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "config/legacy_import.hpp"

namespace {

using namespace ns::config;
namespace fs = std::filesystem;

static const fs::path kInputDir = fs::path(NS_SOURCE_DIR) / "input";

int g_file_counter = 0;

fs::path write_text(const std::string& name, const std::string& text) {
    const auto path = fs::temp_directory_path() /
                      std::format("ns_import_{}_{}.dat", ::getpid(), g_file_counter++);
    std::ofstream out(path);
    out << text;
    return path;
}

TEST(LegacyImport, RampEulerCaseMatchesSourceValues) {
    const auto cfg = import_legacy(kInputDir / "InputRampEuler.cfg");
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;

    EXPECT_EQ(cfg->case_name, "Ramp_KFDS_maxU_24080_2nd");
    EXPECT_EQ(cfg->equations, Equations::Euler);
    EXPECT_EQ(cfg->formulation, Formulation::Nondimensional);
    EXPECT_EQ(cfg->flow, FlowType::Steady);
    EXPECT_EQ(cfg->numerics.time_method, TimeMethod::ForwardEuler);   // time_accuracy = 0
    EXPECT_EQ(cfg->numerics.inviscid_scheme, InviscidScheme::NKFDS_MOVERS_2O);  // 30
    EXPECT_EQ(cfg->numerics.order, SpatialOrder::Second);             // space_accuracy = 1
    EXPECT_DOUBLE_EQ(cfg->numerics.cfl, 0.2);
    EXPECT_DOUBLE_EQ(cfg->grid.scaling, 0.0381);
    EXPECT_EQ(cfg->grid.file, "RampGrid24080.dat");
    EXPECT_EQ(cfg->grid.nx, 240);
    EXPECT_EQ(cfg->grid.ny, 80);
    EXPECT_EQ(cfg->run.max_iterations, 300000);
    EXPECT_DOUBLE_EQ(cfg->physics.mach_inf, 2.0);

    // Topology (1=bottom/jmin, 2=right/imax, 3=top/jmax, 4=left/imin):
    //   10 @4 [2,81]   prescribed-inflow (freestream; no state row)
    //   30 @1 [2,241]  slip-wall
    //   20 @2 [2,81]   transmissive
    //   30 @3 [2,241]  slip-wall
    ASSERT_EQ(cfg->boundaries.size(), 4u);
    EXPECT_TRUE(cfg->boundaries[0].type == BcType::PrescribedInflow);
    EXPECT_FALSE(cfg->boundaries[0].prescribed_state.has_value());  // freestream inflow
    EXPECT_EQ(cfg->boundaries[0].edge, Edge::XMin);
    EXPECT_EQ(cfg->boundaries[0].end, 81);
    EXPECT_TRUE(cfg->boundaries[1].type == BcType::SlipWall);
    EXPECT_EQ(cfg->boundaries[1].edge, Edge::YMin);
    EXPECT_EQ(cfg->boundaries[1].end, 241);
    EXPECT_TRUE(cfg->boundaries[2].type == BcType::Transmissive);
    EXPECT_EQ(cfg->boundaries[2].edge, Edge::XMax);
    EXPECT_TRUE(cfg->boundaries[3].type == BcType::SlipWall);
    EXPECT_EQ(cfg->boundaries[3].edge, Edge::YMax);
}

TEST(LegacyImport, BlasiusNavierStokesCaseMatchesSourceValues) {
    const auto cfg = import_legacy(kInputDir / "InputBlasius.cfg");
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;

    EXPECT_EQ(cfg->case_name, "BlasiusNavMaPt5");
    EXPECT_EQ(cfg->equations, Equations::NavierStokes);
    EXPECT_EQ(cfg->numerics.inviscid_scheme, InviscidScheme::ROE_2O);  // code 24
    EXPECT_EQ(cfg->grid.nx, 160);
    EXPECT_EQ(cfg->grid.ny, 64);
    EXPECT_DOUBLE_EQ(cfg->physics.mach_inf, 0.5);
    EXPECT_DOUBLE_EQ(cfg->physics.re_inf, 5000.0);

    ASSERT_EQ(cfg->boundaries.size(), 5u);
    EXPECT_TRUE(cfg->boundaries[0].type == BcType::Symmetry);            // jmin [2,33]
    EXPECT_EQ(cfg->boundaries[0].start, 2);
    EXPECT_EQ(cfg->boundaries[0].end, 33);
    EXPECT_TRUE(cfg->boundaries[1].type == BcType::Farfield);            // imin
    EXPECT_TRUE(cfg->boundaries[2].type == BcType::Farfield);            // imax
    EXPECT_TRUE(cfg->boundaries[3].type == BcType::Farfield);            // jmax
    EXPECT_TRUE(cfg->boundaries[4].type == BcType::NoSlipAdiabaticWall); // jmin plate
    EXPECT_EQ(cfg->boundaries[4].start, 34);
    EXPECT_EQ(cfg->boundaries[4].end, 161);
}

TEST(LegacyImport, HypersonicHalfCylinderMatchesCommittedConfig) {
    const auto cfg = import_legacy(kInputDir / "InputHypersonicFlowEuler.cfg");
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;

    EXPECT_EQ(cfg->case_name, "HyperCylM20_KfdsMaxU_1st");
    EXPECT_EQ(cfg->formulation, Formulation::Nondimensional);
    EXPECT_EQ(cfg->numerics.time_method, TimeMethod::SSPRK3);          // time_accuracy = 2
    EXPECT_EQ(cfg->numerics.inviscid_scheme, InviscidScheme::NKFDS_MOVERS);  // code 7
    EXPECT_EQ(cfg->numerics.order, SpatialOrder::First);
    EXPECT_EQ(cfg->grid.nx, 80);
    EXPECT_EQ(cfg->grid.ny, 160);
    EXPECT_DOUBLE_EQ(cfg->physics.mach_inf, 20.0);
    ASSERT_TRUE(cfg->physics.reference_state.has_value());
    EXPECT_DOUBLE_EQ(cfg->physics.reference_state->pressure, 100000.0);
}

TEST(LegacyImport, PrescribedStateAndIsothermalWallRows) {
    // Synthetic case exercising the two optional topology row types:
    // an isothermal wall with Twall and a prescribed-inflow with rho,u,v,p.
    const auto dir = fs::temp_directory_path();
    const auto cfg_path = dir / "ns_import_pair_test.cfg";
    {
        std::ofstream out(cfg_path);
        out << "# synthetic legacy pair\n";
        out << "TestIso\n";
        out << "1\n";        // NS
        out << "1\n";        // dimensional -> no reference block
        out << "0\n";        // steady
        out << "0\n";        // forward-euler
        out << "5\n";        // roe-tv
        out << "0\n";        // first order
        out << "0\n";        // green-gauss
        out << "0.4\n";
        out << "1.0\n";
        out << "0\n";        // no restart
        out << "FakeGrid.dat\n";
        out << "20 20\n";
        out << "GridTopFake.dat\n";
        out << "junk_restart.dat\n";
        out << "10\n10\n100\n30\n";
        out << "50000\n0.8\n1.0\n0.0\n101325\n288.15\n";
    }
    {
        std::ofstream out(dir / "GridTopFake.dat");
        out << "# header junk ignored\n# more junk\n";
        out << "2\n";                       // n_seg
        out << "1\n";                       // n_presc
        out << "0\t50\t4\t1\t2\t21\t310.5\n";              // isothermal left, Twall
        out << "1\t10\t2\t21\t2\t21\t1.25\t2.5\t0.0\t9.75\n";  // prescribed right w/ state
    }

    const auto cfg = import_legacy(cfg_path);
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;
    ASSERT_EQ(cfg->boundaries.size(), 2u);

    EXPECT_TRUE(cfg->boundaries[0].type == BcType::NoSlipIsothermalWall);
    EXPECT_EQ(cfg->boundaries[0].edge, Edge::XMin);
    ASSERT_TRUE(cfg->boundaries[0].wall_temperature.has_value());
    EXPECT_DOUBLE_EQ(*cfg->boundaries[0].wall_temperature, 310.5);

    EXPECT_TRUE(cfg->boundaries[1].type == BcType::PrescribedInflow);
    ASSERT_TRUE(cfg->boundaries[1].prescribed_state.has_value());
    EXPECT_DOUBLE_EQ((*cfg->boundaries[1].prescribed_state)[0], 1.25);
    EXPECT_DOUBLE_EQ((*cfg->boundaries[1].prescribed_state)[1], 2.5);
    EXPECT_DOUBLE_EQ((*cfg->boundaries[1].prescribed_state)[3], 9.75);
}

TEST(LegacyImport, ImportedConfigSurvivesTomlRoundTrip) {
    const auto imported = import_legacy(kInputDir / "InputRampEuler.cfg");
    ASSERT_TRUE(imported.has_value());

    const auto reloaded = load(write_text("roundtrip.toml", to_toml(*imported)));
    ASSERT_TRUE(reloaded.has_value()) << "\n" << to_toml(*imported);

    EXPECT_EQ(reloaded->case_name, imported->case_name);
    EXPECT_EQ(reloaded->equations, imported->equations);
    EXPECT_EQ(reloaded->numerics.inviscid_scheme, imported->numerics.inviscid_scheme);
    EXPECT_DOUBLE_EQ(reloaded->grid.scaling, imported->grid.scaling);
    ASSERT_EQ(reloaded->boundaries.size(), imported->boundaries.size());
    for (std::size_t i = 0; i < imported->boundaries.size(); ++i)
        EXPECT_EQ(reloaded->boundaries[i].type, imported->boundaries[i].type);
}

TEST(LegacyImport, MissingFileIsCleanError) {
    const auto cfg =
        import_legacy(fs::temp_directory_path() / "no_such_legacy_cfg.cfg");
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("not found"), std::string::npos);
}

}  // namespace
