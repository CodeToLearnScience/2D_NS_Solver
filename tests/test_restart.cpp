// Restart continuation test: a solver restarted from a mid-run checkpoint
// must reproduce the uninterrupted run's conserved state bit-for-bit
// (serial, deterministic).

#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "config/config.hpp"
#include "mesh/mesh.hpp"
#include "solver/euler_solver.hpp"

namespace {

using namespace ns;

struct CasePaths {
    std::string cfg = std::string(NS_SOURCE_DIR) + "/configs/RampEuler.toml";
    std::string grid = std::string(NS_SOURCE_DIR) + "/input/RampGrid24080.dat";
};

solver::EulerSolver make_solver(config::Config& cfg) {
    auto mesh = StructuredMesh::load_legacy(
        std::string(NS_SOURCE_DIR) + "/input/RampGrid24080.dat", cfg.grid.nx, cfg.grid.ny,
        2, 1.0);
    if (!mesh) throw std::runtime_error(mesh.error().message);
    return solver::EulerSolver(*mesh, cfg);
}

TEST(RestartContinuation, SerialBitwiseContinuation) {
    auto cfg = config::load(CasePaths{}.cfg);
    ASSERT_TRUE(cfg.has_value()) << cfg.error().message;

    const auto tmp = std::filesystem::temp_directory_path();
    const auto rst = tmp / ("ns_cont_" + std::to_string(::getpid()) + ".rst");

    // --- reference run: straight to 16 iterations ---------------------------
    auto ref = make_solver(*cfg);
    ref.run(16);

    // --- split run: 7 iters -> checkpoint -> load -> 9 more -----------------
    auto part = make_solver(*cfg);
    part.run(7);
    part.write_restart(rst);

    auto cont = make_solver(*cfg);          // fresh solver, freestream init
    cont.load_restart(rst);                 // overwrite from checkpoint
    cont.run(9);

    const auto& a = ref.conserved();
    const auto& b = cont.conserved();
    for (int k = 0; k < 4; ++k)
        for (int i = -2; i < a.ni() + 2; ++i)
            for (int j = -2; j < a.nj() + 2; ++j)
                EXPECT_DOUBLE_EQ(a(k, i, j), b(k, i, j))
                    << "plane " << k << " at (" << i << "," << j << ")";

    std::error_code ec;
    std::filesystem::remove(rst, ec);
}

}  // namespace
