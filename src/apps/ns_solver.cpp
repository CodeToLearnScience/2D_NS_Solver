// New-stack driver: runs a steady-Euler forward-Euler case from a canonical
// TOML config and writes legacy-format outputs for parity comparison.
//
//   ns_solver <config.toml> <grid-file> <out-dir> [max-iterations]
//
// Launch under mpirun when built with NS_ENABLE_MPI; every rank then owns a
// j-slab of the domain and rank 0 writes the aggregated outputs.

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "config/config.hpp"
#include "mesh/mesh.hpp"
#include "parallel/decomposition.hpp"
#include "solver/euler_solver.hpp"

#ifdef NS_WITH_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
#ifdef NS_WITH_MPI
    MPI_Init(&argc, &argv);
    int rank = 0, nranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
#else
    const int rank = 0, nranks = 1;
#endif

    int rc = 0;
    try {
        if (argc < 4 || argc > 7) {
            if (rank == 0)
                std::cerr << "usage: ns_solver <config.toml> <grid.dat> <out-dir> "
                             "[max-iters] [--restart <file>]\n";
            return 2;
        }

        const auto cfg = ns::config::load(argv[1]);
        if (!cfg) throw std::runtime_error(cfg.error().message);
        if (cfg->grid.ny < nranks)
            throw std::runtime_error("more ranks than cell rows is unsupported");

        // Full-grid metadata on every rank (file is small; avoids parallel IO
        // for now -- parallel HDF5/XDMF is the documented later upgrade).
        auto full = ns::StructuredMesh::load_legacy(argv[2], cfg->grid.nx,
                                                    cfg->grid.ny, 0, /*scaling=*/1.0);
        if (!full) throw std::runtime_error(full.error().message);

        auto decomp = ns::parallel::SlabDecomposition::make(cfg->grid.ny, nranks, rank);
        auto local = ns::StructuredMesh::load_legacy_slab(
            argv[2], cfg->grid.nx, cfg->grid.ny, decomp.j0, decomp.j1, 2,
            /*scaling=*/1.0);
        if (!local) throw std::runtime_error(local.error().message);

        ns::solver::FullGridMeta meta;
        meta.nx = cfg->grid.nx;
        meta.ny = cfg->grid.ny;
        meta.node_x.assign(static_cast<std::size_t>(meta.nx + 1) * (meta.ny + 1), 0.0);
        meta.node_y.assign(meta.node_x.size(), 0.0);
        for (int j = 0; j <= meta.ny; ++j)
            for (int i = 0; i <= meta.nx; ++i) {
                meta.node_x[static_cast<std::size_t>(j) * (meta.nx + 1) + i] =
                    full->node_x()(i, j);
                meta.node_y[static_cast<std::size_t>(j) * (meta.nx + 1) + i] =
                    full->node_y()(i, j);
            }

        const std::filesystem::path out_dir = argv[3];
        if (rank == 0) std::filesystem::create_directories(out_dir);

        int max_iters = cfg->run.max_iterations;
        std::string restart_file;
        for (int a = 4; a < argc; ++a) {
            const std::string arg = argv[a];
            if (arg == "--restart" && a + 1 < argc) restart_file = argv[++a];
            else max_iters = std::stoi(arg);
        }

        ns::solver::EulerSolver solver(*local, *cfg, decomp, std::move(meta));
        if (!restart_file.empty()) {
            solver.load_restart(restart_file);
            if (rank == 0) std::cout << "restarted from " << restart_file << "\n";
        }
        if (rank == 0) {
            const std::string residue_name = "Residue_" + cfg->case_name +
                                             std::to_string(cfg->grid.nx) + "_" +
                                             std::to_string(cfg->grid.ny) + ".dat";
            solver.open_residue(out_dir / residue_name);
        }

        const int iters = solver.run(max_iters);
        if (rank == 0) std::cout << "completed " << iters << " iterations\n";

        const std::string base = cfg->case_name + "_IsoCont" +
                                 std::to_string(cfg->grid.nx) + "_" +
                                 std::to_string(cfg->grid.ny) + "Iter" +
                                 std::to_string(iters);
        if (std::getenv("NS_DEBUG_PROGRESS"))
            std::fprintf(stderr, "[main r%d writers]\n", rank);
        solver.write_solution(out_dir / (base + ".dat"));
        if (std::getenv("NS_DEBUG_PROGRESS"))
            std::fprintf(stderr, "[main r%d sol-done]\n", rank);
        solver.write_vtk(out_dir / (base + ".vtk"));
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "ns_solver: " << e.what() << "\n";
        rc = 1;
    }

#ifdef NS_WITH_MPI
    MPI_Finalize();
#endif
    return rc;
}
