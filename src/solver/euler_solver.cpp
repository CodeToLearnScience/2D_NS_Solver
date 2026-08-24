#include "solver/euler_solver.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "physics/eos.hpp"

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
      tstep_(mesh.nx(), mesh.ny(), 2),
      ifaces_(kNConv, mesh.nx() + 1, mesh.ny() + 1, 2),
      jfaces_(kNConv, mesh.nx() + 1, mesh.ny() + 1, 2),
      segments_(runtime_segments(cfg.boundaries)) {
    init_freestream();
    // legacy main.cpp: Dependent_Variables BEFORE the pre-loop boundary pass,
    // so the whole canvas starts from cv-reconstructed (not analytic) values
    dependent_variables_all();
    apply_boundaries();
}

void EulerSolver::init_freestream() {
    // Legacy InitFlowNonDimensional, Euler branch.
    rhoinf_ = 1.0;
    tinf_ = 1.0 / kGamma;
    qinf_ = cfg_.physics.mach_inf;
    rgas_ = 1.0;
    cp_ = 1.0 / ((kGamma - 1.0) * cfg_.physics.mach_inf * cfg_.physics.mach_inf);
    pinf_ = rhoinf_ * rgas_ * tinf_;
    alpha_rad_ = cfg_.physics.alpha_deg / (180.0 / (4.0 * std::atan(1.0)));
    uinf_ = qinf_ * std::cos(alpha_rad_);
    vinf_ = qinf_ * std::sin(alpha_rad_);
    ref_visc_ = 0.0;
    eos_.gas_constant = rgas_;  // ND runs run with Rgas = 1 (legacy global)

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

            const double dt = cfg_.numerics.cfl * metrics_.area(i, j) / (sri + srj);
            if (!(dt > 0.0))
                throw std::runtime_error(
                    std::format("non-positive time step at ({},{})", i, j));
            tstep_(i, j) = dt;
        }

    // QUIRK preserved: "local" dt is overridden by the interior minimum.
    double tsmin = 1e32;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (tstep_(i, j) < tsmin) tsmin = tstep_(i, j);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) tstep_(i, j) = tsmin;
}

void EulerSolver::compute_fluxes() {
    diss_.fill(0.0);
    define_convective_faces(mesh_, metrics_, dv_, eos_, ifaces_, jfaces_);

    primitive_differences(dv_, dui_, duj_);

    // QUIRK: cp_ carries Forces()' clobbered value from the previous iteration.
    nkfds_movers_dissipation_second_order(mesh_, metrics_, dv_, dui_, duj_, eos_,
                                          mah_cp_, diss_);

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
void EulerSolver::bc_transmissive(const BoundarySegmentRuntime& s) {
    const int lo = s.start - 2, hi = s.end - 2;
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
            for (int j = lo; j <= hi; ++j)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, nx, j) = dv_(k, nx - 1, j);
                    dv_(k, nx + 1, j) = dv_(k, nx, j);
                }
            break;
        case config::Edge::XMin:
            for (int j = lo; j <= hi; ++j)
                for (int k = 0; k < 8; ++k) {
                    dv_(k, -1, j) = dv_(k, 0, j);
                    dv_(k, -2, j) = dv_(k, -1, j);
                }
            break;
    }
}

void EulerSolver::bc_prescribed_inflow(const BoundarySegmentRuntime& s) {
    double rho = rhoinf_, uu = uinf_, vv = vinf_, pp = pinf_;
    if (s.state) {
        rho = (*s.state)[0];
        uu = (*s.state)[1];
        vv = (*s.state)[2];
        pp = (*s.state)[3];
    }
    const int lo = s.start - 2, hi = s.end - 2;
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
            for (int j = lo; j <= hi; ++j) {
                fill(nx, j);
                for (int k = 0; k < 8; ++k) dv_(k, nx + 1, j) = dv_(k, nx, j);
            }
            break;
        case config::Edge::XMin:
            for (int j = lo; j <= hi; ++j) {
                fill(-1, j);
                for (int k = 0; k < 8; ++k) dv_(k, -2, j) = dv_(k, -1, j);
            }
            break;
    }
}

void EulerSolver::bc_slip_wall(const BoundarySegmentRuntime& s) {
    const Field<Vec2>& si = metrics_.si;
    const Field<Vec2>& sj = metrics_.sj;
    const int lo = s.start - 2, hi = s.end - 2;
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
    } else if (s.edge == config::Edge::XMin) {
        for (int j = lo; j <= hi; ++j) {
            const double ds = std::sqrt(si(0, j).x * si(0, j).x + si(0, j).y * si(0, j).y);
            reflect_fill(-1, j, 0, j, si(0, j).x / ds, si(0, j).y / ds);
            for (int k = 0; k < 8; ++k) dv_(k, -2, j) = dv_(k, -1, j);
        }
    } else {  // XMax
        for (int j = lo; j <= hi; ++j) {
            const double ds = std::sqrt(si(nx, j).x * si(nx, j).x + si(nx, j).y * si(nx, j).y);
            reflect_fill(nx, j, nx - 1, j, -si(nx, j).x / ds, -si(nx, j).y / ds);
            for (int k = 0; k < 8; ++k) dv_(k, nx + 1, j) = dv_(k, nx, j);
        }
    }
}

void EulerSolver::apply_segment(const BoundarySegmentRuntime& s) {
    switch (s.type) {
        case config::BcType::Transmissive: bc_transmissive(s); break;
        case config::BcType::PrescribedInflow: bc_prescribed_inflow(s); break;
        case config::BcType::SlipWall: bc_slip_wall(s); break;
        default: throw std::runtime_error("BC type not supported by EulerSolver yet");
    }
}

void EulerSolver::update_corners() {
    // legacy boundary_conditions.cpp: average the four corner cells' dv across
    // their inner neighbours, then Dependent_Variables_One recomputes them from
    // the (never-updated, still-freestream) cv ghosts.
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
    for (const auto& s : segments_) apply_segment(s);
}

void EulerSolver::scale_rhs_and_update() {
    const int nx = mesh_.nx(), ny = mesh_.ny();
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double adtv = tstep_(i, j) / metrics_.area(i, j);
            for (int k = 0; k < 4; ++k) rhs_(k, i, j) *= adtv;
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
        if (wall->edge == config::Edge::YMin || wall->edge == config::Edge::YMax) {
            const bool bottom = wall->edge == config::Edge::YMin;
            for (int n = wall->start - 2; n <= wall->end - 2; ++n) {
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
            for (int n = wall->start - 2; n <= wall->end - 2; ++n) {
                const int j = n;
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

int EulerSolver::run(int max_iterations) {
    double dt_prev = 0.0;
    int iter = 0;
    while (iter < max_iterations) {
        // QUIRK: legacy accumulates time with the PREVIOUS iteration's dt.
        time_ += dt_prev;
        ++iter;

        cvold_ = cv_;
        diss_.fill(0.0);

        time_step_euler();
        compute_fluxes();
        scale_rhs_and_update();
        dependent_variables_all();
        apply_boundaries();

        const double dt_now = tstep_(0, 0);  // global min after override
        residue_and_forces(iter, dt_now);
        dt_prev = dt_now;

        iter_done_ = iter;
        if (resn1_ < 1e-14) return iter;  // converged
    }
    return iter_done_;
}

void EulerSolver::open_residue(const std::filesystem::path& path) {
    residue_.open(path);
    residue_.flags(std::ios::dec | std::ios::scientific);
    residue_.precision(10);
    if (!residue_)
        throw std::runtime_error("cannot open residue file: " + path.string());
}

void EulerSolver::write_solution(const std::filesystem::path& path) const {
    std::ofstream o(path);
    o.flags(std::ios::dec | std::ios::scientific);
    o.precision(16);
    if (!o) throw std::runtime_error("cannot open solution file: " + path.string());

    const int nx = mesh_.nx(), ny = mesh_.ny();
    o << "TITLE = \"Iter=" << iter_done_ << ", time=" << precise16(time_) << "\"\n"
      << "VARIABLES = \"xc\", \"yc\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"Mach\"\n";
    o << "Zone T=\"" << precise16(time_) << "\", I=" << nx + 1 << ", J=" << ny + 1 << "\n";

    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            auto avg = [&](int k) {
                return 0.25 * (dv_(k, i, j) + dv_(k, i - 1, j) + dv_(k, i - 1, j - 1) +
                               dv_(k, i, j - 1));
            };
            const double rho = avg(RHO), u = avg(U), v = avg(V), p = avg(P),
                         tt = avg(T), a = avg(A);
            const double mach = std::sqrt(u * u + v * v) / a;
            o << mesh_.node_x()(i, j) << "\t" << mesh_.node_y()(i, j) << "\t" << rho
              << "\t" << u << "\t" << v << "\t" << p << "\t" << tt << "\t" << mach
              << "\n";
        }
}

void EulerSolver::write_vtk(const std::filesystem::path& path) const {
    std::ofstream o(path);
    o.flags(std::ios::dec | std::ios::scientific);
    o.precision(16);
    if (!o) throw std::runtime_error("cannot open vtk file: " + path.string());

    const int nx = mesh_.nx(), ny = mesh_.ny();
    const long long n_nodes = static_cast<long long>(nx + 1) * (ny + 1);

    o << "# vtk DataFile Version 3.0\n";
    o << "Field at T=" << precise16(time_) << "\n";
    o << "ASCII\nDATASET STRUCTURED_GRID\n";
    o << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " 1\n";
    o << "POINTS " << n_nodes << " double\n";

    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            o << mesh_.node_x()(i, j) << "\t" << mesh_.node_y()(i, j) << "\t" << 1 << "\n";

    auto cell_avg = [&](int k, int i, int j) {
        return 0.25 * (dv_(k, i, j) + dv_(k, i - 1, j) + dv_(k, i - 1, j - 1) +
                       dv_(k, i, j - 1));
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
