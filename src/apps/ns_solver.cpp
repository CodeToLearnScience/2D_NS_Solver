// New-stack driver: runs a steady-Euler forward-Euler case from a canonical
// TOML config and writes legacy-format outputs for byte-parity comparison.
//
//   ns_solver <config.toml> <grid-file> <out-dir> [max-iterations]

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "config/config.hpp"
#include "mesh/mesh.hpp"
#include "solver/euler_solver.hpp"

int main(int argc, char** argv) {
    if (argc < 4 || argc > 5) {
        std::cerr << "usage: ns_solver <config.toml> <grid.dat> <out-dir> [max-iters]\n";
        return 2;
    }

    try {
        const auto cfg = ns::config::load(argv[1]);
        if (!cfg) throw std::runtime_error(cfg.error().message);

        auto mesh = ns::StructuredMesh::load_legacy(
            argv[2], cfg->grid.nx, cfg->grid.ny, 2, /*scaling=*/1.0);
        if (!mesh) throw std::runtime_error(mesh.error().message);

        const std::filesystem::path out_dir = argv[3];
        std::filesystem::create_directories(out_dir);

        const int max_iters =
            (argc == 5) ? std::stoi(argv[4]) : cfg->run.max_iterations;

        ns::solver::EulerSolver solver(*mesh, *cfg);
        const std::string residue_name = "Residue_" + cfg->case_name +
                                         std::to_string(cfg->grid.nx) + "_" +
                                         std::to_string(cfg->grid.ny) + ".dat";
        solver.open_residue(out_dir / residue_name);

        const int iters = solver.run(max_iters);
        std::cout << "completed " << iters << " iterations\n";

        const std::string base = cfg->case_name + "_IsoCont" +
                                 std::to_string(cfg->grid.nx) + "_" +
                                 std::to_string(cfg->grid.ny) + "Iter" +
                                 std::to_string(iters);
        solver.write_solution(out_dir / (base + ".dat"));
        solver.write_vtk(out_dir / (base + ".vtk"));
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ns_solver: " << e.what() << "\n";
        return 1;
    }
}
