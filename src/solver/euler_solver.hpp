#pragma once

// Single-rank steady-Euler forward-Euler driver (Phase 6).
//
// Reproduces the legacy per-iteration pipeline bit-for-bit:
//   cvold = cv; diss = 0
//   -> Time_Step_Euler (fixed normals; global min-dt override)
//   -> FluxConvDefns -> primitive differences -> NKFDS-MOVERS-2O dissipation
//   -> Avg_Conv_Flux1 RHS assembly
//   -> rhs *= dt/area; cv = cvold - rhs   (+ NaN / negative-pressure guards)
//   -> Dependent_Variables (full canvas)
//   -> Boundary conditions (corner averaging + segment dispatch)
//   -> Residue + forces (including the legacy global-Cp clobber quirk)
//
// Known legacy behaviours intentionally preserved for parity are flagged with
// QUIRK comments throughout.

#include <expected>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "config/config.hpp"
#include "mesh/mesh.hpp"
#include "numerics/face_data.hpp"
#include "numerics/schemes.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"
#include "physics/dependent_variables.hpp"

namespace ns::solver {

struct BoundarySegmentRuntime {
    config::BcType type;
    config::Edge edge;
    int start = 0;  // legacy node index (shifted internally)
    int end = 0;
    std::optional<std::array<double, 4>> state;  // prescribed rho,u,v,p
};

struct FullGridMeta {
    int nx = 0, ny = 0;
    std::vector<double> node_x, node_y;  // (nx+1)*(ny+1), j-major rows
};

class EulerSolver {
public:
    /// Serial / single-patch construction.
    EulerSolver(const StructuredMesh& mesh, const config::Config& cfg);
    /// MPI construction: local patch + decomposition + full-grid metadata.
    EulerSolver(const StructuredMesh& local_mesh, const config::Config& cfg,
                parallel::SlabDecomposition decomp, FullGridMeta meta);

    /// Runs the iteration loop; returns number of iterations performed.
    int run(int max_iterations);

    /// Writes legacy-format outputs into `dir` (residue file is opened once).
    void open_residue(const std::filesystem::path& path);
    void write_solution(const std::filesystem::path& path) const;
    void write_vtk(const std::filesystem::path& path) const;

    [[nodiscard]] const MultiField<double>& conserved() const noexcept { return cv_; }
    [[nodiscard]] const MultiField<double>& primitives() const noexcept { return dv_; }
    // TEMPORARY debug hooks (removed once Phase-6 parity lands)
    [[nodiscard]] const MultiField<double>& dbg_rhs() const noexcept { return rhs_; }
    [[nodiscard]] const MultiField<double>& dbg_diss() const noexcept { return diss_; }
    void dbg_scale_and_update_only() { scale_rhs_and_update(); }
    void dbg_stage_update() { scale_rhs_and_update(); dependent_variables_all(); apply_boundaries(); }
    [[nodiscard]] Field<double> const& dbg_tstep() const noexcept { return tstep_; }
    [[nodiscard]] double dbg_area(int i, int j) const noexcept { return metrics_.area(i, j); }
    void dbg_one_iteration() {
        cvold_ = cv_;
        diss_.fill(0.0);
        time_step_euler();
        compute_fluxes();
    }
    [[nodiscard]] double time() const noexcept { return time_; }

private:
    void init_freestream();
    void fill_meta_from_mesh();
    void finish_init();
    [[nodiscard]] bool is_root() const noexcept { return !mpi_mode_ || decomp_.rank == 0; }
    /// Root-only: assembles a GLOBAL GHOST-FRAMED plane for field k:
    /// size (nx+2)*(ny+2), entry (i,j) at [(i+1)*(ny+2)+(j+1)] covering
    /// i,j in [-1..nx]/[-1..ny]. Exact byte parity with the serial dv_ view
    /// at every sampled node-average location.
    void gather_plane(int k, std::vector<double>& framed) const;
    [[nodiscard]] physics::TransportScaling scaling() const;
    void time_step_euler();
    void compute_fluxes();
    void apply_boundaries();
    void apply_segment(const BoundarySegmentRuntime& seg);
    void update_corners();
    void scale_rhs_and_update();
    void residue_and_forces(int iter, double dt);
    void dependent_variables_all();

    // Boundary kernels (legacy Bc_Transmitive / Bc_Prescribed_Inflow / BC_wall)
    void bc_transmissive(const BoundarySegmentRuntime& s);
    void bc_prescribed_inflow(const BoundarySegmentRuntime& s);
    void bc_slip_wall(const BoundarySegmentRuntime& s);
    void transport_at(int i, int j);
    [[nodiscard]] int decomp_offset() const;
    [[nodiscard]] bool owns_rows(int glo, int ghi) const;

    const StructuredMesh& mesh_;
    const MeshMetrics metrics_;
    const config::Config cfg_;

    MultiField<double> cv_, cvold_, dv_, diss_, rhs_, dui_, duj_;
    Field<double> tstep_;
    numerics::FaceFluxData ifaces_, jfaces_;
    physics::IdealGas eos_{};

    // freestream (ND, Euler branch)
    double rhoinf_ = 0, uinf_ = 0, vinf_ = 0, pinf_ = 0, tinf_ = 0, qinf_ = 0;
    double rgas_ = 1.0, cp_ = 0.0, ref_visc_ = 0.0;
    /// QUIRK: legacy Mahalanobis reads the SEPARATE lowercase global `cp`
    /// that starts at 0 and is clobbered by Forces() every iteration.
    double mah_cp_ = 0.0;
    double alpha_rad_ = 0.0;

    // residue state
    double resn1_ = 0, cl_ = 0, cd_ = 0, cm_ = 0;
    int iresmax_ = 0, jresmax_ = 0;
    double dcv11_ = 0;

    double time_ = 0.0;
    int iter_done_ = 0;
    parallel::SlabDecomposition decomp_{};
    bool mpi_mode_ = false;
    parallel::HaloExchange halo_{};
    FullGridMeta meta_{};
    /// Cached global ghost-framed planes from the last write_solution gather
    /// (avoids a second MPI exchange for write_vtk).
    mutable std::vector<double> cached_planes_[8];
    mutable bool cached_valid_ = false;
    std::ofstream residue_;
    std::vector<BoundarySegmentRuntime> segments_;
};

/// Builds runtime segments from config (legacy node indices preserved).
[[nodiscard]] std::vector<BoundarySegmentRuntime> runtime_segments(
    const std::vector<config::BoundarySegment>& in);

}  // namespace ns::solver
