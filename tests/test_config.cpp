// Unit tests for the TOML config loader and validator (ns::config::load).

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "config/config.hpp"

namespace {

using namespace ns::config;
namespace fs = std::filesystem;

int g_file_counter = 0;

// Writes TOML text to a unique temp file and returns its path.
fs::path write_config(const std::string& text) {
    const auto path = fs::temp_directory_path() /
                      std::format("ns_cfg_test_{}_{}.toml", ::getpid(), g_file_counter++);
    std::ofstream out(path);
    out << text;
    return path;
}

const char* kValidMinimal = R"(# minimal valid steady Euler case
[case]
name = "test-case"
equations = "euler"
formulation = "nondimensional"
flow = "steady"

[grid]
file = "grid.dat"
nx = 10
ny = 10

[numerics]
inviscid_scheme = "llf"
order = "first"
time_method = "ssprk3"
cfl = 0.4

[physics]
Re_inf = 1000.0
Mach_inf = 2.0
p_inf = 101325.0
T_inf = 298.0

[physics.reference_state]
pressure = 100000.0
temperature = 300.0
density = 1.12
velocity = 1.0

[run]
max_iterations = 100

[[boundary]]
edge = "jmin"
start = 2
end = 11
type = "farfield"

[[boundary]]
edge = "imax"
start = 2
end = 11
type = "slip-wall"
)";

TEST(ConfigLoad, MinimalValidConfigParsesWithDefaults) {
    const auto path = write_config(kValidMinimal);
    const auto cfg = load(path);
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;

    EXPECT_EQ(cfg->case_name, "test-case");
    EXPECT_EQ(cfg->equations, Equations::Euler);
    EXPECT_EQ(cfg->formulation, Formulation::Nondimensional);
    EXPECT_EQ(cfg->flow, FlowType::Steady);

    EXPECT_EQ(cfg->grid.file, "grid.dat");
    EXPECT_EQ(cfg->grid.nx, 10);
    EXPECT_EQ(cfg->grid.ny, 10);
    EXPECT_DOUBLE_EQ(cfg->grid.scaling, 1.0);  // default (integer-less file)

    EXPECT_EQ(cfg->numerics.inviscid_scheme, InviscidScheme::LLF);
    EXPECT_EQ(cfg->numerics.order, SpatialOrder::First);
    EXPECT_EQ(cfg->numerics.time_method, TimeMethod::SSPRK3);
    EXPECT_DOUBLE_EQ(cfg->numerics.cfl, 0.4);
    EXPECT_EQ(cfg->numerics.viscous_method, ViscousMethod::GreenGauss);   // default
    EXPECT_DOUBLE_EQ(cfg->numerics.residual_smoothing_eps, 0.0);          // default

    EXPECT_DOUBLE_EQ(cfg->physics.gamma, 1.4);                            // legacy defaults
    EXPECT_DOUBLE_EQ(cfg->physics.prandtl, 0.72);
    EXPECT_DOUBLE_EQ(cfg->physics.gas_constant, 286.92);
    ASSERT_TRUE(cfg->physics.reference_state.has_value());
    EXPECT_DOUBLE_EQ(cfg->physics.reference_state->density, 1.12);

    EXPECT_EQ(cfg->run.max_iterations, 100);
    EXPECT_FALSE(cfg->run.restart);

    EXPECT_EQ(cfg->output.directory, "output");                           // default
    EXPECT_EQ(cfg->output.display_frequency, 100);                        // default
    EXPECT_EQ(cfg->output.write_frequency, 100);                          // default

    ASSERT_EQ(cfg->boundaries.size(), 2u);
    EXPECT_EQ(cfg->boundaries[0].edge, Edge::YMin);
    EXPECT_EQ(cfg->boundaries[0].type, BcType::Farfield);
    EXPECT_EQ(cfg->boundaries[0].start, 2);
    EXPECT_EQ(cfg->boundaries[0].end, 11);
    EXPECT_EQ(cfg->boundaries[1].type, BcType::SlipWall);
}

TEST(ConfigLoad, MissingFileIsCleanError) {
    const auto cfg = load(fs::temp_directory_path() / "ns_cfg_does_not_exist_42.toml");
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("not found"), std::string::npos);
}

TEST(ConfigLoad, MissingRequiredKeyIsReported) {
    std::string text = kValidMinimal;
    text.replace(text.find("cfl = 0.4"), 9, "# cfl removed");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("numerics.cfl"), std::string::npos);
}

TEST(ConfigLoad, UnknownKeyIsRejected) {
    std::string text = kValidMinimal;
    text.replace(text.find("cfl = 0.4"), 9, "cfl = 0.4\nturbo_mode = true");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value()) << "unknown keys must not be silently ignored";
    EXPECT_NE(cfg.error().message.find("numerics.turbo_mode"), std::string::npos);
}

TEST(ConfigLoad, InvalidEnumValueListsAllowedOptions) {
    std::string text = kValidMinimal;
    text.replace(text.find("\"llf\""), 5, "\"warp-drive\"");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("inviscid scheme"), std::string::npos);
    EXPECT_NE(cfg.error().message.find("'llf'"), std::string::npos);
}

TEST(ConfigLoad, MalformedTomlReportsParseError) {
    const auto path = write_config("[case\nname = oops");
    const auto cfg = load(path);
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("TOML parse error"), std::string::npos);
}

TEST(ConfigValidate, SteadyRequiresMaxIterations) {
    std::string text = kValidMinimal;
    text.replace(text.find("max_iterations = 100"), 20, "max_iterations = 0");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("max_iterations"), std::string::npos);
}

TEST(ConfigValidate, UnsteadyRequiresTotalTime) {
    std::string text = kValidMinimal;
    text.replace(text.find("flow = \"steady\""), 15, "flow = \"unsteady\"");
    const auto cfg = load(write_config(text));  // no total_time given
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("total_time"), std::string::npos);
}

TEST(ConfigValidate, RestartWithoutFileRejected) {
    std::string text = kValidMinimal;
    text.replace(text.find("max_iterations = 100"), 20,
                 "max_iterations = 100\nrestart = true");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("restart_file"), std::string::npos);
}

TEST(ConfigValidate, NondimensionalNeedsReferenceState) {
    std::string text = kValidMinimal;
    text.replace(text.find("[physics.reference_state]"),
                 strlen("[physics.reference_state]\npressure = 100000.0\ntemperature = "
                        "300.0\ndensity = 1.12\nvelocity = 1.0\n\n"),
                 "");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("reference_state"), std::string::npos);
}

TEST(ConfigValidate, IsothermalWallNeedsTemperature) {
    std::string text = kValidMinimal;
    text.replace(text.find("type = \"slip-wall\""), 18,
                 "type = \"no-slip-isothermal-wall\"");  // 18 = len("type = \"slip-wall\"")
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("wall_temperature"), std::string::npos);
}

TEST(ConfigValidate, PrescribedInflowWithoutStateIsValid) {
    // Legacy pres_input_flag=0 semantics: no state row -> free stream.
    std::string text = kValidMinimal;
    text.replace(text.find("type = \"farfield\""), 17, "type = \"prescribed-inflow\"");
    const auto cfg = load(write_config(text));
    EXPECT_TRUE(cfg.has_value()) << cfg.error().message;
}

TEST(ConfigValidate, StateOnNonPrescribedSegmentRejected) {
    // Prepend a segment that carries a state but lacks edge/start/end/type;
    // it must be rejected regardless of the valid segments after it.
    std::string text = kValidMinimal;
    text.insert(text.find("[[boundary]]"),
                std::string("[[boundary]]\nstate = [1.0, 2.0, 3.0, 4.0]\n"));
    const auto cfg = load(write_config(text));
    EXPECT_FALSE(cfg.has_value());
}

TEST(ConfigValidate, SegmentEndBeyondEdgeRejected) {
    std::string text = kValidMinimal;                       // nx = ny = 10 -> limit 11
    text.replace(text.find("end = 11"), 8, "end = 99");     // first segment only
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("exceeds the last node"), std::string::npos);
}

TEST(ConfigValidate, MultipleProblemsReportedTogether) {
    std::string text = kValidMinimal;
    text.replace(text.find("cfl = 0.4"), 9, "cfl = -1");
    text.replace(text.find("Mach_inf = 2.0"), 14, "Mach_inf = 0.0");
    const auto cfg = load(write_config(text));
    ASSERT_FALSE(cfg.has_value());
    EXPECT_NE(cfg.error().message.find("cfl"), std::string::npos);
    EXPECT_NE(cfg.error().message.find("Mach_inf"), std::string::npos);
}

TEST(ConfigRoundTrip, ToTomlSurvivesReload) {
    const auto original = load(write_config(kValidMinimal));
    ASSERT_TRUE(original.has_value());

    const std::string text = to_toml(*original);
    const auto reloaded = load(write_config(text));
    ASSERT_TRUE(reloaded.has_value()) << "\n" << text;

    EXPECT_EQ(reloaded->case_name, original->case_name);
    EXPECT_EQ(reloaded->equations, original->equations);
    EXPECT_EQ(reloaded->formulation, original->formulation);
    EXPECT_EQ(reloaded->flow, original->flow);
    EXPECT_EQ(reloaded->numerics.inviscid_scheme, original->numerics.inviscid_scheme);
    EXPECT_EQ(reloaded->numerics.time_method, original->numerics.time_method);
    EXPECT_DOUBLE_EQ(reloaded->numerics.cfl, original->numerics.cfl);
    EXPECT_DOUBLE_EQ(reloaded->grid.scaling, original->grid.scaling);
    ASSERT_EQ(reloaded->boundaries.size(), original->boundaries.size());
    for (std::size_t i = 0; i < original->boundaries.size(); ++i) {
        EXPECT_EQ(reloaded->boundaries[i].edge, original->boundaries[i].edge);
        EXPECT_EQ(reloaded->boundaries[i].type, original->boundaries[i].type);
        EXPECT_EQ(reloaded->boundaries[i].start, original->boundaries[i].start);
        EXPECT_EQ(reloaded->boundaries[i].end, original->boundaries[i].end);
    }
}

}  // namespace
