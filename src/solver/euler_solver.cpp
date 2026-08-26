#include "solver/euler_solver.hpp"

#include <cstdio>
#include <cstdlib>
#include <limits>

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "io/restart.hpp"
#include "physics/gradients.hpp"
#include "physics/viscous_flux.hpp"
#include "physics/eos.hpp"
#include "numerics/llf.hpp"
#include "numerics/schemes.hpp"

namespace ns::solver {
using namespace numerics;
using numerics::kNConv;
namespace {
using physics::cvi::E;
using physics::dvi::A;
using physics::dvi::K;
using physics::dvi::MU;
using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::T;
using physics::dvi::U;
using physics::dvi::V;

constexpr double kGamma = 1.4;  // legacy global Gamma
constexpr double kPr = 0.72;    // legacy global Pr
constexpr double kSuthT0 = 288.15, kSuthS = 110.4;

// Legacy transport fill shared by every BC kernel.
void fill_transport(MultiField<double>& dv, int i, int j, double rgas,
                    physics::Formulation f, double re_inf, double mach_inf,
                    double ref_visc, double cp) {
    dv(T, i, j) = dv(P, i, j) / (dv(RHO, i, j) * rgas);
    dv(A, i, j) = std::sqrt(kGamma * dv(P, i, j) / dv(RHO, i, j));
    if (f == physics::Formulation::Nondimensional) {
        const double c_ratio = kSuthS / kSuthT0;
        dv(MU, i, j) =
            ((1.0 + c_ratio) / (c_ratio + dv(T, i, j))) * std::pow(dv(T, i, j), 1.5) /
            re_inf;
        dv(K, i, j) = dv(MU, i, j) / ((kGamma - 1.0) * mach_inf * mach_inf * kPr);
    } else {
        dv(MU, i, j) = ((kSuthT0 + kSuthS) / (kSuthS + dv(T, i, j))) *
                       std::pow((dv(T, i, j) / kSuthT0), 1.5) * ref_visc;
        dv(K, i, j) = dv(MU, i, j) * cp / kPr;
    }
}

std::string precise16(double v) {
    std::ostringstream out;
    out << std::setprecision(16) << v;
    return out.str();
}

}  // namespace

std::vector<BoundarySegmentRuntime> runtime_segments(
    const std::vector<config::BoundarySegment>& in) {
    std::vector<BoundarySegmentRuntime> out;
    out.reserve(in.size());
    for (const auto& b : in)
        out.push_back({b.type, b.edge, b.start, b.end, b.prescribed_state});
    return out;
}

EulerSolver::EulerSolver(const StructuredMesh& mesh, const config::Config& cfg)
    : mesh_(mesh),
      metrics_(compute_metrics(mesh)),
      cfg_(cfg),
      cv_(4, mesh.nx(), mesh.ny(), 2),
      cvold_(4, mesh.nx(), mesh.ny(), 2),
      dv_(8, mesh.nx(), mesh.ny(), 2),
      diss_(4, mesh.nx(), mesh.ny(), 2),
      rhs_(4, mesh.nx(), mesh.ny(), 2),
      dui_(4, mesh.nx(), mesh.ny(), 2),
      duj_(4, mesh.nx(), mesh.ny(), 2),
      gradfi_(6, mesh.nx() + 1, mesh.ny() + 1, 2),
      gradfj_(6, mesh.nx() + 1, mesh.ny() + 1, 2),
      cui_(4, mesh.nx(), mesh.ny(), 2),
      cuj_(4, mesh.nx(), mesh.ny(), 2),
      tstep_(mesh.nx(), mesh.ny(), 2),
      ifaces_(kNConv, mesh.nx() + 1, mesh.ny() + 1, 2),
      jfaces_(kNConv, mesh.nx() + 1, mesh.ny() + 1, 2),
      segments_(runtime_segments(cfg.boundaries)) {
    fill_meta_from_mesh();
    finish_init();
}

EulerSolver::EulerSolver(const StructuredMesh& local_mesh, const config::Config& cfg,
                         parallel::SlabDecomposition decomp, FullGridMeta meta)
    : mesh_(local_mesh),
      metrics_(compute_metrics(local_mesh)),
      cfg_(cfg),
      cv_(4, local_mesh.nx(), local_mesh.ny(), 2),
      cvold_(4, local_mesh.nx(), local_mesh.ny(), 2),
      dv_(8, local_mesh.nx(), local_mesh.ny(), 2),
      diss_(4, local_mesh.nx(), local_mesh.ny(), 2),
      rhs_(4, local_mesh.nx(), local_mesh.ny(), 2),
      dui_(4, local_mesh.nx(), local_mesh.ny(), 2),
      duj_(4, local_mesh.nx(), local_mesh.ny(), 2),
      gradfi_(6, local_mesh.nx() + 1, local_mesh.ny() + 1, 2),
      gradfj_(6, local_mesh.nx() + 1, local_mesh.ny() + 1, 2),
      cui_(4, local_mesh.nx(), local_mesh.ny(), 2),
      cuj_(4, local_mesh.nx(), local_mesh.ny(), 2),
      tstep_(local_mesh.nx(), local_mesh.ny(), 2),
      ifaces_(kNConv, local_mesh.nx() + 1, local_mesh.ny() + 1, 2),
      jfaces_(kNConv, local_mesh.nx() + 1, local_mesh.ny() + 1, 2),
      segments_(runtime_segments(cfg.boundaries)),
      decomp_(decomp),
      mpi_mode_(decomp.nranks > 1),
      halo_(local_mesh.nx(), decomp.ny_local(), 2,
            decomp.owns_lower_neighbor() ? decomp.rank - 1 : -1,
            decomp.owns_upper_neighbor() ? decomp.rank + 1 : -1),
      meta_(std::move(meta)) {
    finish_init();
}

void EulerSolver::setup_mpi_gradient_tables() {
    const int nx = mesh_.nx(), ny = mesh_.ny();

    // 1) Corrected sj face metrics: copy, then exchange ghost face rows so
    //    slots -1 and ny+1 hold the TRUE neighbour face vectors.
    sj_corrected_ = metrics_.sj;
#ifdef NS_WITH_MPI
    halo_.exchange_face_metric_rows(&sj_corrected_.data()->x, &sj_corrected_.data()->y,
                                    sj_corrected_.stride_i() * 2, ny + 1);
#endif

    // 2) True neighbour cell areas via a one-time plane exchange.
    Field<double> area_x(metrics_.area);
#ifdef NS_WITH_MPI
    halo_.exchange(area_x);
#endif

    // 3) Per-row j-face normalization with two-sided interface dual volumes.
    rvol_j_ = Field<double>(nx, ny + 1, 2);
    const bool bot_global = halo_.at_global_bottom();
    const bool top_global = halo_.at_global_top();
    for (int i = 0; i < nx; ++i) {
        for (int jf = 1; jf < ny; ++jf)
            rvol_j_(i, jf) = 2.0 / (metrics_.area(i, jf) + metrics_.area(i, jf - 1));
        rvol_j_(i, 0) = bot_global ? 2.0 / metrics_.area(i, 0)
                                   : 2.0 / (area_x(i, -1) + metrics_.area(i, 0));
        rvol_j_(i, ny) = top_global ? 2.0 / metrics_.area(i, ny - 1)
                                    : 2.0 / (metrics_.area(i, ny - 1) + area_x(i, ny));
    }
}

void EulerSolver::fill_meta_from_mesh() {
    meta_.nx = mesh_.nx();
    meta_.ny = mesh_.ny();
    meta_.node_x.assign(static_cast<std::size_t>(meta_.nx + 1) * (meta_.ny + 1), 0.0);
    meta_.node_y.assign(meta_.node_x.size(), 0.0);
    for (int j = 0; j <= meta_.ny; ++j)
        for (int i = 0; i <= meta_.nx; ++i) {
            meta_.node_x[static_cast<std::size_t>(j) * (meta_.nx + 1) + i] =
                mesh_.node_x()(i, j);
            meta_.node_y[static_cast<std::size_t>(j) * (meta_.nx + 1) + i] =
                mesh_.node_y()(i, j);
        }
}

void EulerSolver::finish_init() {
    init_freestream();
#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        halo_.exchange(dv_);
        setup_mpi_gradient_tables();
    }
#endif
    // legacy main.cpp: Dependent_Variables BEFORE the pre-loop boundary pass,
    // so the whole canvas starts from cv-reconstructed (not analytic) values
    dependent_variables_all();
#ifdef NS_WITH_MPI
    if (mpi_mode_) halo_.exchange(dv_);
#endif
    apply_boundaries();
}

void EulerSolver::init_freestream() {
    if (cfg_.formulation == config::Formulation::Dimensional) {
        // Legacy InitFlowDimensional.
        rgas_ = cfg_.physics.gas_constant > 0 ? cfg_.physics.gas_constant : 287.0;
        rhoinf_ = cfg_.physics.p_inf / (rgas_ * cfg_.physics.t_inf);
        qinf_ = cfg_.physics.mach_inf *
                std::sqrt(kGamma * rgas_ * cfg_.physics.t_inf);
        cp_ = kGamma * rgas_ / (kGamma - 1.0);
        pinf_ = cfg_.physics.p_inf;
        tinf_ = cfg_.physics.t_inf;
        ref_visc_ = cfg_.physics.re_inf > 0
                        ? 1.458e-6 * std::pow(tinf_, 1.5) / (tinf_ + kSuthS)
                        : 0.0;
        alpha_rad_ = cfg_.physics.alpha_deg / (180.0 / (4.0 * std::atan(1.0)));
        uinf_ = qinf_ * std::cos(alpha_rad_);
        vinf_ = qinf_ * std::sin(alpha_rad_);
        eos_.gas_constant = rgas_;

        const int ng = cv_.ng();
        const double c_ratio = kSuthS / kSuthT0;
        for (int i = -ng; i < cv_.ni() + ng; ++i)
            for (int j = -ng; j < cv_.nj() + ng; ++j) {
                dv_(RHO,i,j)=rhoinf_; dv_(U,i,j)=uinf_; dv_(V,i,j)=vinf_;
                dv_(P,i,j)=pinf_; dv_(T,i,j)=tinf_;
                dv_(A,i,j)=std::sqrt(kGamma*pinf_/rhoinf_);
                dv_(MU,i,j)=((kSuthT0+kSuthS)/(kSuthS+tinf_))*
                             std::pow(tinf_/kSuthT0,1.5)*ref_visc_;
                dv_(K,i,j)=dv_(MU,i,j)*cp_/kPr;
                cv_(RHO,i,j)=rhoinf_;
                cv_(U,i,j)=rhoinf_*uinf_;
                cv_(V,i,j)=rhoinf_*vinf_;
                cv_(E,i,j)=pinf_/0.4+0.5*rhoinf_*(uinf_*uinf_+vinf_*vinf_);
            }
        return;
    }

    // Nondimensional (legacy InitFlowNonDimensional).  Legacy splits the
    // scaling by equation type -- Euler: T*=1/gamma, q*=M, R*=1;
    // Navier-Stokes: T*=q*=1, R*=1/(gamma*M^2).  Prescribed-inflow states
    // (e.g. SWBLI post-shock data) are expressed in these internal units, so
    // the split is mandatory for unit consistency.
    rhoinf_ = 1.0;
    cp_ = 1.0 / ((kGamma - 1.0) * cfg_.physics.mach_inf * cfg_.physics.mach_inf);
    alpha_rad_ = cfg_.physics.alpha_deg / (180.0 / (4.0 * std::atan(1.0)));
    const double mach2 = cfg_.physics.mach_inf * cfg_.physics.mach_inf;
    if (cfg_.equations == config::Equations::NavierStokes) {
        tinf_ = 1.0;
        qinf_ = 1.0;
        rgas_ = 1.0 / (kGamma * mach2);
        ref_visc_ = cfg_.physics.reference_length > 0
                        ? cfg_.physics.reference_length / cfg_.physics.re_inf
                        : 1.0 / cfg_.physics.re_inf;
    } else {
        tinf_ = 1.0 / kGamma;
        qinf_ = cfg_.physics.mach_inf;
        rgas_ = 1.0;
        ref_visc_ = 0.0;
    }
    pinf_ = rhoinf_ * rgas_ * tinf_;
    uinf_ = qinf_ * std::cos(alpha_rad_);
    vinf_ = qinf_ * std::sin(alpha_rad_);
    eos_.gas_constant = rgas_;

    const int ng = cv_.ng();
    for (int i = -ng; i < cv_.ni() + ng; ++i)
        for (int j = -ng; j < cv_.nj() + ng; ++j) {
            dv_(RHO, i, j) = rhoinf_;
            dv_(U, i, j) = uinf_;
            dv_(V, i, j) = vinf_;
            dv_(P, i, j) = pinf_;
            dv_(T, i, j) = tinf_;
            dv_(A, i, j) = std::sqrt(kGamma * pinf_ / rhoinf_);
            const double c_ratio = kSuthS / kSuthT0;
            dv_(MU, i, j) = ((1.0 + c_ratio) / (c_ratio + tinf_)) *
                            std::pow(tinf_, 1.5) / cfg_.physics.re_inf;
            dv_(K, i, j) = dv_(MU, i, j) / ((kGamma - 1.0) * cfg_.physics.mach_inf *
                                            cfg_.physics.mach_inf * kPr);

            cv_(RHO, i, j) = rhoinf_;
            cv_(U, i, j) = rhoinf_ * uinf_;
            cv_(V, i, j) = rhoinf_ * vinf_;
            cv_(E, i, j) =
                pinf_ / (kGamma - 1.0) + 0.5 * rhoinf_ * (uinf_ * uinf_ + vinf_ * vinf_);
        }
}

physics::TransportScaling EulerSolver::scaling() const {
    physics::TransportScaling s{};
    s.re_inf = cfg_.physics.re_inf;
    s.mach_inf = cfg_.physics.mach_inf;
    s.ref_visc = ref_visc_;
    s.cp = cp_;
    return s;
}

void EulerSolver::time_step_euler() {
    const int nx = mesh_.nx(), ny = mesh_.ny();
    const Field<Vec2>& si = metrics_.si;
    const Field<Vec2>& sj = metrics_.sj;

    tstep_.fill(0.0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double u = dv_(U, i, j), v = dv_(V, i, j);

            const double sxi = 0.5 * (si(i, j).x + si(i + 1, j).x);
            const double syi = 0.5 * (si(i, j).y + si(i + 1, j).y);
            const double dsi = std::sqrt(sxi * sxi + syi * syi);
            const double vci = u * (sxi / dsi) + v * (syi / dsi);
            const double sri = (std::fabs(vci) + dv_(A, i, j)) * dsi;

            const double sxj = 0.5 * (sj(i, j).x + sj(i, j + 1).x);
            const double syj = 0.5 * (sj(i, j).y + sj(i, j + 1).y);
            const double dsj = std::sqrt(sxj * sxj + syj * syj);
            const double vcj = u * (sxj / dsj) + v * (syj / dsj);
            const double srj = (std::fabs(vcj) + dv_(A, i, j)) * dsj;

            double sr_total = sri + srj;

            // Viscous spectral radius (legacy Time_Step_NS): adds diffusion
            // limit dt ∝ Δx²·ρ/μ. Only for Navier-Stokes.
#ifdef NS_WITH_MPI
            const bool is_ns = cfg_.equations == config::Equations::NavierStokes;
            if (!mpi_mode_ || is_ns)
#endif
            if (cfg_.equations == config::Equations::NavierStokes) {
                constexpr double cfac = 2.0;
                const double rho_inv = 1.0 / dv_(RHO, i, j);
                const double f1 = (4.0 / 3.0) * rho_inv;
                const double f2 = kGamma * rho_inv;
                const double fac = std::max(f1, f2);
                const double mu = dv_(MU, i, j);
                const double dtv = fac * mu / eos_.prandtl;
                const double ar_ = metrics_.area(i, j);
                const double srvi = cfac * dtv * dsi * dsi / ar_;
                const double srvj = cfac * dtv * dsj * dsj / ar_;
                sr_total += srvi + srvj;
            }

            const double dt = cfg_.numerics.cfl * metrics_.area(i, j) / sr_total;
            if (!(dt > 0.0))
                throw std::runtime_error(
                    std::format("non-positive time step at ({},{})", i, j));
            tstep_(i, j) = dt;
#ifdef NS_WITH_MPI
#endif
            if (std::getenv("NS_DEBUG_DT") && i==0 && j==0)
                std::fprintf(stderr, "[DT] sri=%.6g srj=%.6g sr_total=%.6g area=%.6g dt=%.6g mu=%.6g\n",
                             sri, srj, sr_total, metrics_.area(i,j), dt, dv_(MU,i,j));
        }

    // QUIRK preserved: "local" dt is overridden by the interior minimum.
    // MPI: min is order-independent, so dt stays bitwise-identical to serial.
    double tsmin = 1e32;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (tstep_(i, j) < tsmin) tsmin = tstep_(i, j);
#ifdef NS_WITH_MPI
    if (mpi_mode_)
        MPI_Allreduce(MPI_IN_PLACE, &tsmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) tstep_(i, j) = tsmin;
}

void EulerSolver::compute_fluxes() {
    diss_.fill(0.0);
    define_convective_faces(mesh_, metrics_, dv_, eos_, ifaces_, jfaces_);

    primitive_differences(dv_, dui_, duj_);

    // QUIRK: mah_cp_ carries Forces()' clobbered value from the previous iteration.
    double emax_hint = std::numeric_limits<double>::quiet_NaN();
    const bool is_second_order = cfg_.numerics.order == config::SpatialOrder::Second;

    switch (cfg_.numerics.inviscid_scheme) {
        case config::InviscidScheme::NKFDS_MOVERS_2O:
#ifdef NS_WITH_MPI
            if (mpi_mode_) {
                emax_hint = numerics::compute_emax_window(mesh_, dv_, eos_);
                MPI_Allreduce(MPI_IN_PLACE, &emax_hint, 1, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
            }
#endif
            nkfds_movers_dissipation_second_order(mesh_, metrics_, dv_, dui_, duj_,
                                                  eos_, mah_cp_, diss_, emax_hint);
            break;
        case config::InviscidScheme::NKFDS_MOVERS:
            nkfds_movers_dissipation_first_order(mesh_, metrics_, dv_, eos_, mah_cp_,
                                                 diss_);
            break;
        case config::InviscidScheme::ROE_2O:
            roe_dissipation_second_order(mesh_, metrics_, dv_, dui_, duj_, eos_, diss_);
            break;
        case config::InviscidScheme::LLF:
            numerics::llf_dissipation_first_order(mesh_, metrics_, dv_, ifaces_.con_var_diff,
                                        jfaces_.con_var_diff, eos_, diss_);
            break;
        case config::InviscidScheme::LLF_2O:
            numerics::llf_dissipation_second_order(mesh_, metrics_, dv_, dui_, duj_,
                                                   eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS:
            numerics::movers_dissipation_muscl(mesh_, metrics_, dv_, dui_, duj_,
                                               eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS_H1:
            numerics::movers_h1_dissipation(mesh_, metrics_, dv_, eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS_LE1:
            conserved_differences(cv_, cui_, cuj_);
            numerics::movers_le1_dissipation(mesh_, metrics_, cv_, dv_, cui_, cuj_,
                                             eos_, diss_);
            break;
        case config::InviscidScheme::ROE_TV:
            numerics::roe_tv_dissipation(mesh_, metrics_, dv_, eos_, diss_);
            break;
        case config::InviscidScheme::KFDS:
            kfds_dissipation_first_order(mesh_, metrics_, dv_, ifaces_, jfaces_,
                                         eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS_H2:
            movers_h2_dissipation(mesh_, metrics_, dv_, dui_, duj_, 0.5, 0.5,
                                  eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS_LE2:
            conserved_differences(cv_, cui_, cuj_);
            numerics::movers_le2_dissipation(mesh_, metrics_, cv_, dv_, cui_, cuj_,
                                             dui_, duj_, eos_, diss_);
            break;
        case config::InviscidScheme::ROE_TV_2O:
            numerics::roe_tv_2o_dissipation(mesh_, metrics_, dv_, dui_, duj_,
                                            eos_, diss_);
            break;
        case config::InviscidScheme::KFDS_2O:
            kfds_dissipation_second_order(mesh_, metrics_, dv_, ifaces_, jfaces_,
                                          eos_, diss_);
            break;
        case config::InviscidScheme::ECCS:
            eccs_dissipation(mesh_, metrics_, dv_, dui_, duj_, eos_, diss_);
            break;
        case config::InviscidScheme::MOVERS_NWSC:
            movers_nwsc_dissipation(mesh_, metrics_, dv_, dui_, duj_, eos_, diss_);
            break;
        default:
            throw std::runtime_error("scheme not yet wired into EulerSolver");
    }

    // Viscous fluxes for Navier-Stokes
    if (cfg_.equations == config::Equations::NavierStokes) {
        if (mpi_mode_) {
            const physics::FaceGradBoundaryFlags flg{halo_.at_global_bottom(),
                                                    halo_.at_global_top()};
            const physics::FaceGradTables tbl{&sj_corrected_, &rvol_j_};
            physics::green_gauss_face_gradients(mesh_, metrics_, dv_, gradfi_,
                                                gradfj_, flg, &tbl);
        } else {
            physics::green_gauss_face_gradients(mesh_, metrics_, dv_, gradfi_,
                                                gradfj_);
        }
        physics::accumulate_viscous_flux(mesh_, metrics_, dv_, gradfi_, gradfj_, diss_);
    }

    assemble_rhs_from_average_fluxes(mesh_, metrics_, ifaces_.avg_flux,
                                     jfaces_.avg_flux, diss_, rhs_);
}

void EulerSolver::transport_at(int i, int j) {
    fill_transport(dv_, i, j, rgas_, physics::Formulation::Nondimensional,
                   cfg_.physics.re_inf, cfg_.physics.mach_inf, ref_visc_, cp_);
}

// ---------------------------------------------------------------------------
// Boundary kernels -- legacy index translation: bnode/sbind/ebind are legacy
// node indices; subtract 2 for the new canvas.
// ---------------------------------------------------------------------------
int EulerSolver::decomp_offset() const {
#ifdef NS_WITH_MPI
    return mpi_mode_ ? decomp_.j0 : 0;
#else
    return 0;
#endif
}
bool EulerSolver::owns_rows(int glo, int ghi) const {
#ifdef NS_WITH_MPI
    if (!mpi_mode_) return true;
    return std::max(glo, decomp_.j0) <= std::min(ghi, decomp_.j1 - 1);
#else
    (void)glo; (void)ghi;
    return true;
#endif
}

void EulerSolver::bc_transmissive(const BoundarySegmentRuntime& s) {
#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        if (s.edge == config::Edge::YMin && !decomp_.owns_jmin()) return;
        if (s.edge == config::Edge::YMax && !decomp_.owns_jmax()) return;
    }
#endif
    int lo = s.start - 2, hi = s.end - 2;
    const int nx = mesh_.nx(), ny = mesh_.ny();

    switch (s.edge) {
        case config::Edge::YMin:
            for (int i = lo; i <= hi; ++i)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, i, -1) = dv_(k, i, 0);
                    dv_(k, i, -2) = dv_(k, i, -1);
                }
            break;
        case config::Edge::YMax:
            for (int i = lo; i <= hi; ++i)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, i, ny) = dv_(k, i, ny - 1);
                    dv_(k, i, ny + 1) = dv_(k, i, ny);
                }
            break;
        case config::Edge::XMax:
#ifdef NS_WITH_MPI
            lo = mpi_mode_ ? std::max(lo, decomp_.j0) : lo;
            hi = mpi_mode_ ? std::min(hi, decomp_.j1 - 1) : hi;
#endif
            for (int j = lo - decomp_offset(); j <= hi - decomp_offset(); ++j)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, nx, j) = dv_(k, nx - 1, j);
                    dv_(k, nx + 1, j) = dv_(k, nx, j);
                }
            break;
        case config::Edge::XMin:
#ifdef NS_WITH_MPI
            lo = mpi_mode_ ? std::max(lo, decomp_.j0) : lo;
            hi = mpi_mode_ ? std::min(hi, decomp_.j1 - 1) : hi;
#endif
            for (int j = lo - decomp_offset(); j <= hi - decomp_offset(); ++j)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, -1, j) = dv_(k, 0, j);
                    dv_(k, -2, j) = dv_(k, -1, j);
                }
            break;
    }
}

void EulerSolver::bc_prescribed_inflow(const BoundarySegmentRuntime& s) {
#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        if (s.edge == config::Edge::YMin && !decomp_.owns_jmin()) return;
        if (s.edge == config::Edge::YMax && !decomp_.owns_jmax()) return;
    }
#endif
    double rho = rhoinf_, uu = uinf_, vv = vinf_, pp = pinf_;
    if (s.state) {
        rho = (*s.state)[0];
        uu = (*s.state)[1];
        vv = (*s.state)[2];
        pp = (*s.state)[3];
    }
    int lo = s.start - 2, hi = s.end - 2;
    const int nx = mesh_.nx(), ny = mesh_.ny();

    auto fill = [&](int ci, int cj) {
        dv_(RHO, ci, cj) = rho;
        dv_(U, ci, cj) = uu;
        dv_(V, ci, cj) = vv;
        dv_(P, ci, cj) = pp;
        transport_at(ci, cj);
    };

    switch (s.edge) {
        case config::Edge::YMin:
            for (int i = lo; i <= hi; ++i) {
                fill(i, -1);
                for (int k = 0; k < 8; ++k) dv_(k, i, -2) = dv_(k, i, -1);
            }
            break;
        case config::Edge::YMax:
            for (int i = lo; i <= hi; ++i) {
                fill(i, ny);
                for (int k = 0; k < 8; ++k) dv_(k, i, ny + 1) = dv_(k, i, ny);
            }
            break;
        case config::Edge::XMax:
#ifdef NS_WITH_MPI
            lo = mpi_mode_ ? std::max(lo, decomp_.j0) : lo;
            hi = mpi_mode_ ? std::min(hi, decomp_.j1 - 1) : hi;
#endif
            for (int j = lo - decomp_offset(); j <= hi - decomp_offset(); ++j) {
                fill(nx, j);
                for (int k = 0; k < 8; ++k) dv_(k, nx + 1, j) = dv_(k, nx, j);
            }
            break;
        case config::Edge::XMin:
#ifdef NS_WITH_MPI
            lo = mpi_mode_ ? std::max(lo, decomp_.j0) : lo;
            hi = mpi_mode_ ? std::min(hi, decomp_.j1 - 1) : hi;
#endif
            for (int j = lo - decomp_offset(); j <= hi - decomp_offset(); ++j) {
                fill(-1, j);
                for (int k = 0; k < 8; ++k) dv_(k, -2, j) = dv_(k, -1, j);
            }
            break;
    }
}

void EulerSolver::bc_slip_wall(const BoundarySegmentRuntime& s) {
#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        if (s.edge == config::Edge::YMin && !decomp_.owns_jmin()) return;
        if (s.edge == config::Edge::YMax && !decomp_.owns_jmax()) return;
    }
#endif
    const Field<Vec2>& si = metrics_.si;
    const Field<Vec2>& sj = metrics_.sj;
    int lo = s.start - 2, hi = s.end - 2;
    const int nx = mesh_.nx(), ny = mesh_.ny();

    auto reflect_fill = [&](int gi, int gj, int ii, int ij, double nx, double ny) {
        const double vnor0 = dv_(U, ii, ij) * nx + dv_(V, ii, ij) * ny;
        const double vtan = dv_(U, ii, ij) * ny - dv_(V, ii, ij) * nx;
        const double vnor = -vnor0;
        dv_(RHO, gi, gj) = dv_(RHO, ii, ij);
        dv_(U, gi, gj) = vnor * nx + vtan * ny;
        dv_(V, gi, gj) = vnor * ny - vtan * nx;
        dv_(P, gi, gj) = dv_(P, ii, ij);
        transport_at(gi, gj);
    };

    if (s.edge == config::Edge::YMin) {
        for (int i = lo; i <= hi; ++i) {
            const double ds = std::sqrt(sj(i, 0).x * sj(i, 0).x + sj(i, 0).y * sj(i, 0).y);
            reflect_fill(i, -1, i, 0, sj(i, 0).x / ds, sj(i, 0).y / ds);
            for (int k = 0; k < 8; ++k) dv_(k, i, -2) = dv_(k, i, -1);
        }
    } else if (s.edge == config::Edge::YMax) {
        for (int i = lo; i <= hi; ++i) {
            const double ds = std::sqrt(sj(i, ny).x * sj(i, ny).x + sj(i, ny).y * sj(i, ny).y);
            reflect_fill(i, ny, i, ny - 1, -sj(i, ny).x / ds, -sj(i, ny).y / ds);
            for (int k = 0; k < 8; ++k) dv_(k, i, ny + 1) = dv_(k, i, ny);
        }
    } else {  // XMin / XMax
#ifdef NS_WITH_MPI
        lo = mpi_mode_ ? std::max(lo, decomp_.j0) : lo;
        hi = mpi_mode_ ? std::min(hi, decomp_.j1 - 1) : hi;
#endif
        const bool left = s.edge == config::Edge::XMin;
        for (int j = lo - decomp_offset(); j <= hi - decomp_offset(); ++j) {
            if (left) {
                const double ds =
                    std::sqrt(si(0, j).x * si(0, j).x + si(0, j).y * si(0, j).y);
                reflect_fill(-1, j, 0, j, si(0, j).x / ds, si(0, j).y / ds);
                for (int k = 0; k < 8; ++k) dv_(k, -2, j) = dv_(k, -1, j);
            } else {
                const double ds =
                    std::sqrt(si(nx, j).x * si(nx, j).x + si(nx, j).y * si(nx, j).y);
                reflect_fill(nx, j, nx - 1, j, -si(nx, j).x / ds,
                             -si(nx, j).y / ds);
                for (int k = 0; k < 8; ++k) dv_(k, nx + 1, j) = dv_(k, nx, j);
            }
        }
    }
}

void EulerSolver::bc_no_slip_wall(const BoundarySegmentRuntime& s,
                                  std::optional<double> wall_temp) {
    // Isothermal walls override T at ghost with the specified wall temperature
    const bool isothermal = wall_temp.has_value() && *wall_temp > 0.0;
    const Field<Vec2>& si = metrics_.si;
    const Field<Vec2>& sj = metrics_.sj;
    const int lo = s.start - 2, hi = s.end - 2;
    const int nx = mesh_.nx(), ny = mesh_.ny();

    auto fill_wall = [&](int gi, int gj, int ii, int ij) {
        // No-slip: u=v=0 at ghost (velocity reflected with zero magnitude)
        dv_(RHO, gi, gj) = dv_(RHO, ii, ij);
        dv_(U, gi, gj) = -dv_(U, ii, ij);   // reflected normal component
        dv_(V, gi, gj) = -dv_(V, ii, ij);
        dv_(P, gi, gj) = dv_(P, ii, ij);
        transport_at(gi, gj);
    };

    if (s.edge == config::Edge::YMin) {
        for (int i = lo; i <= hi; ++i) {
            fill_wall(i,-1,i,0);
            for (int k = 0; k < 8; ++k) dv_(k,i,-2)=dv_(k,i,-1);
        }
    } else if (s.edge == config::Edge::YMax) {
        for (int i = lo; i <= hi; ++i) {
            fill_wall(i,ny,i,ny-1);
            for (int k = 0; k < 8; ++k) dv_(k,i,ny+1)=dv_(k,i,ny);
        }
    } else if (s.edge == config::Edge::XMin) {
        for (int j = lo; j <= hi; ++j) {
            fill_wall(-1,j,0,j);
            for (int k = 0; k < 8; ++k) dv_(k,-2,j)=dv_(k,-1,j);
        }
    } else {  // XMax
        for (int j = lo; j <= hi; ++j) {
            fill_wall(nx,j,nx-1,j);
            for (int k = 0; k < 8; ++k) dv_(k,nx+1,j)=dv_(k,nx,j);
        }
    }
}

void EulerSolver::bc_farfield(const BoundarySegmentRuntime& s) {
    // Characteristic far-field: use freestream state (simplified but adequate
    // for supersonic external flows where all characteristics enter).
    bc_prescribed_inflow(s);  // same as prescribed inflow with freestream state
}

void EulerSolver::bc_symmetry(const BoundarySegmentRuntime& s) {
#ifdef NS_WITH_MPI
    if (mpi_mode_ && (s.edge == config::Edge::YMin && !decomp_.owns_jmin())) return;
    if (mpi_mode_ && (s.edge == config::Edge::YMax && !decomp_.owns_jmax())) return;
#endif
    const int lo = s.start - 2, hi = s.end - 2;
    const int nx = mesh_.nx(), ny = mesh_.ny();

    auto sym_fill = [&](int gi, int gj, int ii, int ij) {
        for (int k = 0; k < 8; ++k) dv_(k, gi, gj) = dv_(k, ii, ij);
        dv_(V, gi, gj) = -dv_(V, ii, ij);  // negate normal velocity only
    };

    if (s.edge == config::Edge::YMin) {
        for (int i = lo; i <= hi; ++i) {
            sym_fill(i,-1,i,0);
            for (int k = 0; k < 8; ++k) dv_(k,i,-2)=dv_(k,i,-1);
        }
    } else if (s.edge == config::Edge::YMax) {
        for (int i = lo; i <= hi; ++i) {
            sym_fill(i,ny,i,ny-1);
            for (int k = 0; k < 8; ++k) dv_(k,i,ny+1)=dv_(k,i,ny);
        }
    } else if (s.edge == config::Edge::XMin) {
        for (int j = lo; j <= hi; ++j) {
            sym_fill(-1,j,0,j);
            for (int k = 0; k < 8; ++k) dv_(k,-2,j)=dv_(k,-1,j);
        }
    } else {  // XMax
        for (int j = lo; j <= hi; ++j) {
            sym_fill(nx,j,nx-1,j);
            for (int k = 0; k < 8; ++k) dv_(k,nx+1,j)=dv_(k,nx,j);
        }
    }
}

void EulerSolver::clamp_segment_bounds(BoundarySegmentRuntime& s) {
    const int ng = 2;  // cv_.ng()
    const int max_i = mesh_.nx() + ng - 1;
    const int max_j = mesh_.ny() + ng - 1;
    int lo = std::max(s.start - 2, -ng);
    int hi = std::min(s.end - 2,
                      (s.edge == config::Edge::XMin || s.edge == config::Edge::XMax)
                          ? mesh_.ny() + ng - 1 : max_i);
    s.start = lo + 2;
    s.end = hi + 2;
}

void EulerSolver::apply_segment(const BoundarySegmentRuntime& s) {
    switch (s.type) {
        case config::BcType::Transmissive: bc_transmissive(s); break;
        case config::BcType::PrescribedInflow: bc_prescribed_inflow(s); break;
        case config::BcType::SlipWall: bc_slip_wall(s); break;
        case config::BcType::NoSlipAdiabaticWall: bc_no_slip_wall(s); break;
        case config::BcType::NoSlipIsothermalWall: bc_no_slip_wall(s); break;
        case config::BcType::Farfield: bc_farfield(s); break;
        case config::BcType::Cut: bc_transmissive(s); break;
        case config::BcType::Symmetry: bc_symmetry(s); break;
    }
}

void EulerSolver::update_corners() {
    // legacy boundary_conditions.cpp: average the four corner cells' dv across
    // their inner neighbours, then Dependent_Variables_One recomputes them from
    // the (never-updated, still-freestream) cv ghosts.
#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        const bool top = decomp_.owns_jmax(), bottom = decomp_.owns_jmin();
        const int ib = mesh_.nx() - 1, jb = mesh_.ny() - 1;
        auto do_corner = [&](int ci, int cj) {
            const int di = ci < 0 ? 1 : -1;
            const int dj = cj < 0 ? 1 : -1;
            for (int k = 0; k < 8; ++k)
                dv_(k, ci, cj) = 0.5 * (dv_(k, ci + di, cj) + dv_(k, ci, cj + dj));
            compute_dependent_variables_at(cv_, dv_, ci, cj, eos_,
                                           physics::Formulation::Nondimensional,
                                           scaling());
        };
        if (bottom) {
            do_corner(-1, -1);
            do_corner(ib + 1, -1);
        }
        if (top) {
            do_corner(-1, jb + 1);
            do_corner(ib + 1, jb + 1);
        }
        return;
    }
#endif
    const int ib = mesh_.nx() - 1, jb = mesh_.ny() - 1;
    const int corners[4][2] = {{-1, -1}, {ib + 1, -1}, {-1, jb + 1}, {ib + 1, jb + 1}};
    for (auto& c : corners) {
        const int ci = c[0], cj = c[1];
        const int di = ci < 0 ? 1 : -1;
        const int dj = cj < 0 ? 1 : -1;
        for (int k = 0; k < 8; ++k)
            dv_(k, ci, cj) = 0.5 * (dv_(k, ci + di, cj) + dv_(k, ci, cj + dj));
        compute_dependent_variables_at(cv_, dv_, ci, cj, eos_,
                                       physics::Formulation::Nondimensional, scaling());
    }
}

void EulerSolver::apply_boundaries() {
    update_corners();
    for (const auto& s : segments_) {
        // j-slab ownership: only the rank owning the global bottom/top may
        // stamp YMin/YMax segments; otherwise inter-rank ghost rows would be
        // clobbered with physical-BC state instead of neighbour data.
        if (s.edge == config::Edge::YMin && !decomp_.owns_jmin()) continue;
        if (s.edge == config::Edge::YMax && !decomp_.owns_jmax()) continue;
        apply_segment(s);
    }
}

void EulerSolver::scale_rhs_and_update() {
    const int nx = mesh_.nx(), ny = mesh_.ny();
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double adtv = tstep_(i, j) / metrics_.area(i, j);
            for (int k = 0; k < 4; ++k) rhs_(k, i, j) *= adtv;
        }

    if (std::getenv("NS_DEBUG_PROGRESS")) {
        double pmin = 1e300;
        int pi_ = 0, pj_ = 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const double pr = eos_.pressure(cv_(RHO, i, j), cv_(U, i, j),
                                                cv_(V, i, j), cv_(E, i, j));
                if (pr < pmin) { pmin = pr; pi_ = i; pj_ = j; }
            }
        std::fprintf(stderr, "[R%d pre-update minp=%.6g at (%d,%d)] cv@(0,0): %.6g %.6g %.6g %.6g\n",
                     decomp_.rank, pmin, pi_, pj_,
                     cv_(RHO,0,0), cv_(U,0,0), cv_(V,0,0), cv_(E,0,0));
    }
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < 4; ++k)
                cv_(k, i, j) = cvold_(k, i, j) - rhs_(k, i, j);
            const double p =
                (cv_(E, i, j) -
                 0.5 * (cv_(U, i, j) * cv_(U, i, j) + cv_(V, i, j) * cv_(V, i, j)) /
                     cv_(RHO, i, j)) *
                (eos_.gamma - 1.0);
            bool bad = !(p >= 0.0);
            for (int k = 0; k < 4 && !bad; ++k) bad = !(cv_(k, i, j) == cv_(k, i, j));
            if (bad)
                throw std::runtime_error(
                    std::format("NaN/negative state at cell ({},{})", i, j));
        }
}

void EulerSolver::dependent_variables_all() {
    compute_dependent_variables(cv_, dv_, eos_, physics::Formulation::Nondimensional,
                                scaling());
}

void EulerSolver::residue_and_forces(int iter, double dt) {
    const int nx = mesh_.nx(), ny = mesh_.ny();

    // ---- residual norms ---------------------------------------------------
    double rn1 = 0, rn2 = 0, rn3 = 0, rn4 = 0, dcv1max = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double d1 = cv_(RHO, i, j) - cvold_(RHO, i, j);
            rn1 += d1 * d1;
            if (std::fabs(d1) >= dcv1max) {
                dcv1max = std::fabs(d1);
                iresmax_ = i + 2;  // legacy index reported
                jresmax_ = j + 2;
            }
                rn2 += std::pow(cv_(U, i, j) - cvold_(U, i, j), 2);
            rn3 += std::pow(cv_(V, i, j) - cvold_(V, i, j), 2);
            rn4 += std::pow(cv_(E, i, j) - cvold_(E, i, j), 2);
        }

#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        double partial[4] = {rn1, rn2, rn3, rn4};
        MPI_Allreduce(MPI_IN_PLACE, partial, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        rn1 = partial[0];
        rn2 = partial[1];
        rn3 = partial[2];
        rn4 = partial[3];
        // location of |drho| maximum: largest (globalJ,globalI) among ties,
        // approximating the legacy last-wins scan
        double lv = dcv1max;
        int loc = (jresmax_ << 16) | iresmax_;
        MPI_Allreduce(MPI_IN_PLACE, &lv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (lv > 0 && std::fabs(dcv1max) == lv) {
            MPI_Allreduce(MPI_IN_PLACE, &loc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            jresmax_ = loc >> 16;
            iresmax_ = loc & 0xFFFF;
        } else if (lv > 0) {
            MPI_Allreduce(MPI_IN_PLACE, &loc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            jresmax_ = loc >> 16;
            iresmax_ = loc & 0xFFFF;
        }
    }
#endif

    if (iter == 1) {
        dcv11_ = std::sqrt(rn1) + 1e-32;
        resn1_ = 1.0;
    } else {
        resn1_ = std::sqrt(rn1);
    }
    (void)dcv11_;

    // ---- forces (QUIRK: only the LAST wall segment is integrated; global cp
    // is clobbered with the last local pressure coefficient) -----------------
    double cx = 0, cy = 0;
    cm_ = 0;
    const double xref = 0.25, yref = 0.0, cref = 1.0;

    const BoundarySegmentRuntime* wall = nullptr;
    for (const auto& s : segments_)
        if (s.type == config::BcType::SlipWall ||
            s.type == config::BcType::NoSlipAdiabaticWall ||
            s.type == config::BcType::NoSlipIsothermalWall)
            wall = &s;

    if (wall != nullptr) {
        const Field<Vec2>& si = metrics_.si;
        const Field<Vec2>& sj = metrics_.sj;
#ifdef NS_WITH_MPI
        // vertical walls span all j: every rank integrates its slab overlap
        auto clip_lo = [&](int n0) {
            return mpi_mode_ ? std::max(n0, decomp_.j0) : n0;
        };
        auto clip_hi = [&](int n1) {
            return mpi_mode_ ? std::min(n1, decomp_.j1 - 1) : n1;
        };
#endif
        if (wall->edge == config::Edge::YMin || wall->edge == config::Edge::YMax) {
            const bool bottom = wall->edge == config::Edge::YMin;
            bool owner = true;
#ifdef NS_WITH_MPI
            owner = !mpi_mode_ || (bottom ? decomp_.owns_jmin() : decomp_.owns_jmax());
#endif
            for (int n = (owner ? wall->start - 2 : 1); owner && n <= wall->end - 2; ++n) {
                const int i = n;
                const int inner = bottom ? 0 : ny - 1;
                const int gface = bottom ? 0 : ny;
                double sx, sy;
                if (bottom) {
                    sx = sj(i, gface).x;
                    sy = sj(i, gface).y;
                } else {
                    sx = -sj(i, gface).x;
                    sy = -sj(i, gface).y;
                }
                const double pwall = 0.5 * (dv_(P, i, bottom ? -1 : ny) + dv_(P, i, inner));
                const double cp_loc =
                    2.0 * (pwall - pinf_) / (rhoinf_ * qinf_ * qinf_);
                // legacy uses node row/column ins1 (the INTERIOR side), not
                // the ghost-face row
                const int nrow = bottom ? 0 : ny - 1;
                const double xa = (0.5 * (mesh_.node_x()(i, nrow) +
                                          mesh_.node_x()(i + 1, nrow)) - xref) / cref;
                const double ya = (0.5 * (mesh_.node_y()(i, nrow) +
                                          mesh_.node_y()(i + 1, nrow)) - yref) / cref;
                const double dcy = sy * cp_loc;
                const double dcx = sx * cp_loc;
                cy += dcy;
                cx += dcx;
                cm_ += dcx * ya - dcy * xa;
                mah_cp_ = cp_loc;  // QUIRK: clobbered every call
            }
        } else {
            const bool left = wall->edge == config::Edge::XMin;
            const int wlo = wall->start - 2, whi = wall->end - 2;
#ifdef NS_WITH_MPI
            const int wlo2 = mpi_mode_ ? std::max(wlo, decomp_.j0) : wlo;
            const int whi2 = mpi_mode_ ? std::min(whi, decomp_.j1 - 1) : whi;
#else
            const int wlo2 = wlo, whi2 = whi;
#endif
                for (int n = wlo2; n <= whi2; ++n) {
                const int j = n - decomp_offset();  // global -> local row
                const int inner = left ? 0 : nx - 1;
                const int gface = left ? 0 : nx;
                double sx, sy;
                if (left) {
                    sx = si(gface, j).x;
                    sy = si(gface, j).y;
                } else {
                    sx = -si(gface, j).x;
                    sy = -si(gface, j).y;
                }
                const double pwall =
                    0.5 * (dv_(P, left ? -1 : nx, j) + dv_(P, inner, j));
                const double cp_loc = 2.0 * (pwall - pinf_) / (rhoinf_ * qinf_ * qinf_);
                const int ncol = left ? 0 : nx - 1;
                const double xa = (0.5 * (mesh_.node_x()(ncol, j) +
                                          mesh_.node_x()(ncol, j + 1)) - xref) / cref;
                const double ya = (0.5 * (mesh_.node_y()(ncol, j) +
                                          mesh_.node_y()(ncol, j + 1)) - yref) / cref;
                const double dcy = sy * cp_loc;
                const double dcx = sx * cp_loc;
                cy += dcy;
                cx += dcx;
                cm_ += dcx * ya - dcy * xa;
                mah_cp_ = cp_loc;
            }
        }
    }

    cl_ = cy * std::cos(alpha_rad_) - cx * std::sin(alpha_rad_);
    cd_ = cy * std::sin(alpha_rad_) - cx * std::cos(alpha_rad_);

    // ---- residue file ------------------------------------------------------
    if (residue_.is_open()) {
        if (iter == 1)
            residue_ << "TITLE = Residue and forces\n"
                     << "VARIABLES = Iter, dt, tot_t, Res_rho, iresmax, jresmax, Cl, Cd, "
                        "Cm\n";
        residue_ << iter << "\t" << dt << "\t" << time_ << "\t" << resn1_ << "\t"
                 << iresmax_ << "\t" << jresmax_ << "\t" << cl_ << "\t" << cd_ << "\t"
                 << cm_ << "\n";
        residue_.flush();
    }
}

void EulerSolver::scale_rhs_only() {
    const int nx = mesh_.nx(), ny = mesh_.ny();
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double adtv = tstep_(i, j) / metrics_.area(i, j);
            for (int k = 0; k < 4; ++k) rhs_(k, i, j) *= adtv;
        }
}

void EulerSolver::ssprk_blend(double ark, double brk) {
    const int nx = mesh_.nx(), ny = mesh_.ny();
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < 4; ++k) {
                if (ark == 0.0)
                    cv_(k,i,j) = cvold_(k,i,j) - rhs_(k,i,j);
                else
                    cv_(k,i,j) = ark*cvold_(k,i,j) + brk*(cv_(k,i,j)-rhs_(k,i,j));
            }
            const double p = eos_.pressure(
                cv_(RHO,i,j), cv_(U,i,j), cv_(V,i,j), cv_(E,i,j));
            bool bad = !(p >= 0.0);
            for (int k = 0; k < 4 && !bad; ++k) bad = !(cv_(k,i,j)==cv_(k,i,j));
            if (bad)
                throw std::runtime_error(std::format(
                    "NaN/negative state at cell ({},{})", i, j));
        }
}

int EulerSolver::run(int max_iterations) {
    int nstages = 1;
    double ark[3] = {0.0, 0.0, 0.0};
    double brk[3] = {1.0, 1.0, 1.0};
    switch (cfg_.numerics.time_method) {
        case config::TimeMethod::SSPRK2:
            nstages = 2; ark[0]=0;ark[1]=0.5; brk[0]=1;brk[1]=0.5; break;
        case config::TimeMethod::SSPRK3:
            nstages = 3; ark[0]=0;ark[1]=0.75;ark[2]=1.0/3.0;
            brk[0]=1;brk[1]=0.25;brk[2]=2.0/3.0; break;
        default: break;
    }

    double dt_prev = 0.0;
    int iter = 0;
    while (iter < max_iterations) {
        if (cfg_.flow == config::FlowType::Unsteady &&
            time_ >= cfg_.run.total_time) break;

        time_ += dt_prev;
        ++iter;

        cvold_ = cv_;

        time_step_euler();

        for (int s = 0; s < nstages; ++s) {
            diss_.fill(0.0);
            compute_fluxes();
            scale_rhs_only();
            ssprk_blend(ark[s], brk[s]);
            dependent_variables_all();
#ifdef NS_WITH_MPI
            if (mpi_mode_) halo_.exchange(dv_);
#endif
            apply_boundaries();
        }

        const double dt_now = tstep_(0, 0);
        residue_and_forces(iter, dt_now);
        dt_prev = dt_now;
        iter_done_ = iter;
        if (resn1_ < 1e-14) return iter;
    }
    return iter_done_;
}

void EulerSolver::gather_plane(int k, std::vector<double>& framed) const {
    gather_plane_from(dv_, k, framed);
}

void EulerSolver::gather_plane_from(const MultiField<double>& src, int k,
                                    std::vector<double>& framed) const {
    using parallel::SlabDecomposition;
    const bool gp_dbg = std::getenv("NS_DEBUG_PROGRESS") != nullptr;
    if (gp_dbg)
        std::fprintf(stderr, "[G%d r%d enter]\n", k,
#ifdef NS_WITH_MPI
                     decomp_.rank
#else
                     0
#endif
        );
    const int nx = meta_.nx, ny = meta_.ny;
    const int W = ny + 2;
    framed.assign(static_cast<std::size_t>(nx + 2) * W, 0.0);
    if (cached_valid_ && !cached_planes_[k].empty()) {
        framed = cached_planes_[k];
        return;
    }

#ifdef NS_WITH_MPI
    if (mpi_mode_) {
        auto pack = [&](std::vector<double>& buf, int nrows) {
            buf.resize(static_cast<std::size_t>(mesh_.nx() + 4) * (nrows + 4));
            for (int a = 0; a < mesh_.nx() + 4; ++a)
                for (int b = 0; b < nrows + 4; ++b)
                    buf[static_cast<std::size_t>(a) * (nrows + 4) + b] =
                        src(k, a - 2, b - 2);
        };
        auto place = [&](int j0, int nrows, const std::vector<double>& buf) {
            for (int a = 0; a < mesh_.nx() + 4; ++a)
                for (int b = 0; b < nrows + 4; ++b) {
                    const int gi = a - 2, gj = j0 + b - 2;
                    if (gi < -1 || gi > nx || gj < -1 || gj > ny) continue;
                    framed[static_cast<std::size_t>(gi + 1) * W + (gj + 1)] =
                        buf[static_cast<std::size_t>(a) * (nrows + 4) + b];
                }
        };

        std::vector<double> mine;
        pack(mine, decomp_.ny_local());
        if (!is_root()) {
            MPI_Send(mine.data(), static_cast<int>(mine.size()), MPI_DOUBLE, 0,
                     900 + k, MPI_COMM_WORLD);
            if (gp_dbg) std::fprintf(stderr, "[G%d sent]\n", k);
            return;
        }
        if (gp_dbg) std::fprintf(stderr, "[G%d own-placed]\n", k);
        place(decomp_.j0, decomp_.ny_local(), mine);  // root's own slab
        for (int r = 1; r < decomp_.nranks; ++r) {
            auto d = SlabDecomposition::make(decomp_.ny_global, decomp_.nranks, r);
            std::vector<double> buf(static_cast<std::size_t>(mesh_.nx() + 4) *
                                    (d.ny_local() + 4));
            MPI_Recv(buf.data(), static_cast<int>(buf.size()), MPI_DOUBLE, r, 900 + k,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            place(d.j0, d.ny_local(), buf);
        }
        cached_planes_[k] = framed;
        cached_valid_ = true;
        return;
    }
#endif
    if (gp_dbg) std::fprintf(stderr, "[G%d done-serial]\n", k);
    // serial: direct copy of the ghosted view restricted to [-1..nx][-1..ny]
    for (int i = -1; i <= nx; ++i)
        for (int j = -1; j <= ny; ++j)
            framed[static_cast<std::size_t>(i + 1) * W + (j + 1)] = src(k, i, j);
    cached_planes_[k] = framed;
    cached_valid_ = true;
}


void EulerSolver::write_restart(const std::filesystem::path& path) const {
    std::vector<std::vector<double>> framed(4);
    for (int k = 0; k < 4; ++k) gather_plane_from(cv_, k, framed[k]);

    if (!is_root()) return;  // only root owns the file

    if (std::getenv("NS_DEBUG_PROGRESS"))
        std::fprintf(stderr, "[WR root=%d]\n", static_cast<int>(is_root()));
    static constexpr const char* kNames[4] = {"rho", "rhou", "rhov", "E"};
    if (std::getenv("NS_DEBUG_PROGRESS"))
        for (int k = 0; k < 4; ++k)
            std::fprintf(stderr, "[WR k%d min=%.6g max=%.6g]\n", k,
                         *std::min_element(framed[k].begin(), framed[k].end()),
                         *std::max_element(framed[k].begin(), framed[k].end()));
    std::vector<std::pair<std::string, std::vector<double>>> pairs;
    for (int k = 0; k < 4; ++k) pairs.emplace_back(kNames[k], std::move(framed[k]));

    // framed planes are (nx+2)x(ny+2): declare zero extra ghosts so the
    // reader's canvas arithmetic matches byte-for-byte.
    auto r = io::write_restart_values(path, meta_.nx + 2, meta_.ny + 2, 0, pairs);
    if (!r) throw std::runtime_error(r.error().message);
}

void EulerSolver::load_restart(const std::filesystem::path& path) {
#ifdef NS_WITH_MPI
    const bool mpi = mpi_mode_;
#else
    constexpr bool mpi = false;
#endif
    io::RestartData data;
    std::vector<double> sendbuf;

    if (!mpi || is_root()) {
        auto r = io::read_restart(path);
        if (!r) throw std::runtime_error(r.error().message);
        if (r->nx != mesh_.nx() + 2 || r->ny != meta_.ny + 2 || r->ng != 0)
            throw std::runtime_error(std::format(
                "restart dims {}x{} do not match case frame {}x{}", r->nx, r->ny,
                mesh_.nx() + 2, meta_.ny + 2));
        data = std::move(*r);
        // strip ghost frames into a plane-major contiguous real-cell buffer:
        // [plane k][i][j]
        sendbuf.reserve(static_cast<std::size_t>(4) * mesh_.nx() * meta_.ny);
        for (auto& f : data.fields)
            for (int i = 0; i < mesh_.nx(); ++i)
                for (int j = 0; j < meta_.ny; ++j)
                    sendbuf.push_back(f.values[static_cast<std::size_t>(i + 1) *
                                                   (meta_.ny + 2) +
                                               (j + 1)]);
    }

    const std::size_t plane_size =
        static_cast<std::size_t>(mesh_.nx()) * mesh_.ny();

    for (int k = 0; k < 4; ++k) {
        std::vector<double> local(plane_size);
#ifdef NS_WITH_MPI
        if (mpi) {
            std::vector<int> counts(decomp_.nranks), displs(decomp_.nranks);
            int off = 0;
            for (int r = 0; r < decomp_.nranks; ++r) {
                auto d = parallel::SlabDecomposition::make(decomp_.ny_global,
                                                           decomp_.nranks, r);
                counts[r] = mesh_.nx() * d.ny_local();
                displs[r] = off;
                off += counts[r];
            }
            // sendbuf planes are laid out over the GLOBAL row count
            const std::size_t plane_off =
                static_cast<std::size_t>(k) * mesh_.nx() * meta_.ny;
            MPI_Scatterv(is_root() ? sendbuf.data() + plane_off : nullptr,
                         counts.data(), displs.data(), MPI_DOUBLE, local.data(),
                         static_cast<int>(local.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (std::getenv("NS_DEBUG_PROGRESS"))
            std::fprintf(stderr, "[SCAT r%d k%d first=%g %g]\n",
                         decomp_.rank, k, local[0], local[1]);
        } else
#endif
        {
            const io::RestartField& f = data.fields[static_cast<std::size_t>(k)];
            for (int i = 0; i < mesh_.nx(); ++i)
                for (int j = 0; j < meta_.ny; ++j)
                    local[static_cast<std::size_t>(i) * mesh_.ny() + j] =
                        f.values[static_cast<std::size_t>(i + 1) * (meta_.ny + 2) +
                                 (j + 1)];
        }
        for (int i = 0; i < mesh_.nx(); ++i)
            for (int j = 0; j < mesh_.ny(); ++j)
                cv_(k, i, j) = local[static_cast<std::size_t>(i) * mesh_.ny() + j];
    }

    if (std::getenv("NS_DEBUG_PROGRESS")) {
        double pmin = 1e300;
        int pi_ = 0, pj_ = 0;
        for (int j = 0; j < mesh_.ny(); ++j)
            for (int i = 0; i < mesh_.nx(); ++i) {
                const double pr = eos_.pressure(cv_(RHO, i, j), cv_(U, i, j),
                                                cv_(V, i, j), cv_(E, i, j));
                if (pr < pmin) { pmin = pr; pi_ = i; pj_ = j; }
            }
        std::fprintf(stderr, "[LOAD r%d minp=%.6g at (%d,%d)] cv00=%.3g %.3g %.3g %.3g\n",
                     decomp_.rank, pmin, pi_, pj_,
                     cv_(RHO,0,0), cv_(U,0,0), cv_(V,0,0), cv_(E,0,0));
    }

    // mirror finish_init tail: reconstruct everything downstream of cv
    dependent_variables_all();
#ifdef NS_WITH_MPI
    if (mpi_mode_) halo_.exchange(dv_);
#endif
    apply_boundaries();
}

void EulerSolver::open_residue(const std::filesystem::path& path) {
    if (!is_root()) return;  // rank 0 owns the residue file
    residue_.open(path);
    residue_.flags(std::ios::dec | std::ios::scientific);
    residue_.precision(10);
    if (!residue_)
        throw std::runtime_error("cannot open residue file: " + path.string());
}

void EulerSolver::write_solution(const std::filesystem::path& path) const {
    const bool dbgw = std::getenv("NS_DEBUG_PROGRESS") != nullptr;
    if (dbgw)
        std::fprintf(stderr, "[ws enter root=%d]\n", static_cast<int>(is_root()));
    std::vector<double> gp[8];
    for (int k = 0; k < 8; ++k) {
        if (dbgw) std::fprintf(stderr, "[W sol g%d]\n", k);
        gather_plane(k, gp[k]);
        if (dbgw)
            std::fprintf(stderr, "[W sol G%d size=%zu root=%d]\n", k, gp[k].size(),
                         static_cast<int>(is_root()));
    }
    const int nx = meta_.nx, ny = meta_.ny;
    const int W = ny + 2;

    if (!is_root()) return;

    std::ofstream o(path);
    o.flags(std::ios::dec | std::ios::scientific);
    o.precision(16);
    if (!o) throw std::runtime_error("cannot open solution file: " + path.string());

    o << "TITLE = \"Iter=" << iter_done_ << ", time=" << precise16(time_) << "\"\n"
      << "VARIABLES = \"xc\", \"yc\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"Mach\"\n";
    o << "Zone T=\"" << precise16(time_) << "\", I=" << nx + 1 << ", J=" << ny + 1 << "\n";

    auto node_at = [&](char which, int i, int j) -> double {
        const auto& arr = which == 'x' ? meta_.node_x : meta_.node_y;
        return arr[static_cast<std::size_t>(j) * (nx + 1) + i];
    };
    auto cell_at = [&](int k, int i, int j) -> double {
        return gp[k][static_cast<std::size_t>(i + 1) * W + (j + 1)];
    };
    if (dbgw) std::fprintf(stderr, "[ws writing]\n");
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            auto avg = [&](int k) {
                return 0.25 * (cell_at(k, i, j) + cell_at(k, i - 1, j) +
                               cell_at(k, i - 1, j - 1) + cell_at(k, i, j - 1));
            };
            const double rho = avg(RHO), u = avg(U), v = avg(V), p = avg(P),
                         tt = avg(T), a = avg(A);
            const double mach = std::sqrt(u * u + v * v) / a;
            o << node_at('x', i, j) << "\t" << node_at('y', i, j) << "\t" << rho
              << "\t" << u << "\t" << v << "\t" << p << "\t" << tt << "\t" << mach
              << "\n";
        }
}

void EulerSolver::write_vtk(const std::filesystem::path& path) const {
    // Reuses planes gathered by write_solution -- no second exchange, so
    // non-root ranks must not enter at all (nothing to match on the wire).
    if (!is_root() || !cached_valid_) return;
    const int nx = meta_.nx, ny = meta_.ny;
    const int W = ny + 2;
    auto cell_at = [&](int k, int i, int j) -> double {
        return cached_planes_[k][static_cast<std::size_t>(i + 1) * W + (j + 1)];
    };

    std::ofstream o(path);
    o.flags(std::ios::dec | std::ios::scientific);
    o.precision(16);
    if (!o) throw std::runtime_error("cannot open vtk file: " + path.string());
    const long long n_nodes = static_cast<long long>(nx + 1) * (ny + 1);

    o << "# vtk DataFile Version 3.0\n";
    o << "Field at T=" << precise16(time_) << "\n";
    o << "ASCII\nDATASET STRUCTURED_GRID\n";
    o << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " 1\n";
    o << "POINTS " << n_nodes << " double\n";

    auto node_at = [&](char which, int i, int j) -> double {
        const auto& arr = which == 'x' ? meta_.node_x : meta_.node_y;
        return arr[static_cast<std::size_t>(j) * (nx + 1) + i];
    };
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            o << node_at('x', i, j) << "\t" << node_at('y', i, j) << "\t" << 1 << "\n";

    auto cell_avg = [&](int k, int i, int j) {
        return 0.25 * (cell_at(k, i, j) + cell_at(k, i - 1, j) +
                       cell_at(k, i - 1, j - 1) + cell_at(k, i, j - 1));
    };

    o << "POINT_DATA " << n_nodes << "\n";
    o << "SCALARS density double\nLOOKUP_TABLE default\n";
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) o << cell_avg(RHO, i, j) << "\n";

    o << "SCALARS pressure double\nLOOKUP_TABLE default\n";
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) o << cell_avg(P, i, j) << "\n";

    o << "VECTORS velocity double\n";
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            o << cell_avg(U, i, j) << " " << cell_avg(V, i, j) << " " << 0.0 << "\n";
}

}  // namespace ns::solver
