#pragma once

// Typed configuration model for the refactored solver (Phase 2).
//
// A config is loaded from TOML, validated up-front, and passed around as an
// immutable value -- there are no globals anywhere downstream of this header.
//
// Legacy mappings preserved verbatim (see REFACTORING_PLAN.md):
//   eqn_type        -> Equations          (0 euler, 1 navier-stokes)
//   dimen           -> Formulation        (0 nondimensional, 1 dimensional)
//   flow_type       -> FlowType           (0 steady, 1 unsteady)
//   time_accuracy   -> TimeMethod         (0 forward-euler, 1 ssprk2, 2 ssprk3)
//   space_accuracy  -> SpatialOrder       (0 first, 1 second)
//   inviscid_scheme -> InviscidScheme     (see kSchemeNames below)
//   bc flags        -> BcType             (10 prescribed-inflow, 20 transmissive,
//                                          30 slip-wall, 40 no-slip-adiabatic,
//                                          50 no-slip-isothermal, 60 farfield,
//                                          70 symmetry, 80 cut)
//   bound_ind       -> Edge               (1 bottom->jmin, 2 right->imax,
//                                          3 top->jmax, 4 left->imin)

#include <array>
#include <expected>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace ns::config {

enum class Equations { Euler, NavierStokes };
enum class Formulation { Nondimensional, Dimensional };
enum class FlowType { Steady, Unsteady };
enum class TimeMethod { ForwardEuler, SSPRK2, SSPRK3 };
enum class SpatialOrder { First, Second };
enum class ViscousMethod { GreenGauss, LeastSquares };

enum class Edge { XMin, XMax, YMin, YMax };

enum class BcType {
    PrescribedInflow,        // legacy flag 10
    Transmissive,            // 20
    SlipWall,                // 30 (inviscid / Euler wall)
    NoSlipAdiabaticWall,     // 40
    NoSlipIsothermalWall,    // 50
    Farfield,                // 60
    Symmetry,                // 70
    Cut                      // 80
};

enum class InviscidScheme {
    LLF,              // legacy 1st-order code 0
    MOVERS,           // 1  (also invoked by 2nd-order code 28)
    MOVERS_H1,        // 2
    MOVERS_LE1,       // 3
    ROE_TV,           // 5
    KFDS,             // 6
    NKFDS_MOVERS,     // 7
    LLF_2O,           // 2nd-order code 20
    MOVERS_H2,        // 22
    MOVERS_LE2,       // 23
    ROE_2O,           // 24
    ROE_TV_2O,        // 25
    KFDS_2O,          // 26
    ECCS,             // 27
    MOVERS_NWSC,      // 29
    NKFDS_MOVERS_2O   // 30
};

struct BoundarySegment {
    Edge edge{};
    int start = 0;    // first node index along the edge (legacy convention: interior nodes start at 2)
    int end = 0;      // last node index (inclusive)
    BcType type{};
    std::optional<std::array<double, 4>> prescribed_state;  // {rho, u, v, p}
    std::optional<double> wall_temperature;                 // isothermal walls only
};

struct GridConfig {
    std::string file;
    int nx = 0;
    int ny = 0;
    double scaling = 1.0;
};

struct NumericsConfig {
    InviscidScheme inviscid_scheme{};
    SpatialOrder order{};
    TimeMethod time_method = TimeMethod::ForwardEuler;
    double cfl = 0.0;
    ViscousMethod viscous_method = ViscousMethod::GreenGauss;
    double residual_smoothing_eps = 0.0;  // legacy epsirs; unused until Phase 6
};

struct PhysicsConfig {
    double gamma = 1.4;
    double prandtl = 0.72;
    double gas_constant = 286.92;   // legacy Rgas
    double re_inf = 0.0;
    double mach_inf = 0.0;
    double alpha_deg = 0.0;
    double p_inf = 0.0;
    double t_inf = 0.0;
    double reference_length = 1.0;  // legacy Lref

    struct ReferenceState {
        double pressure = 0.0;
        double temperature = 0.0;
        double density = 0.0;
        double velocity = 0.0;
    };
    std::optional<ReferenceState> reference_state;  // nondimensional runs only
};

struct RunControl {
    int max_iterations = 0;
    double total_time = 0.0;      // stop time for unsteady runs
    bool restart = false;
    std::string restart_file;     // resolved against the output directory (Phase 6)
};

struct OutputConfig {
    std::string directory = "output";
    int display_frequency = 100;
    int write_frequency = 100;
    // valid entries: "solution", "vtk", "surface", "residual"
    std::vector<std::string> formats = {"solution", "vtk", "surface", "residual"};
};

struct Config {
    std::string case_name;
    Equations equations{};
    Formulation formulation{};
    FlowType flow{};
    GridConfig grid;
    NumericsConfig numerics;
    PhysicsConfig physics;
    RunControl run;
    OutputConfig output;
    std::vector<BoundarySegment> boundaries;
};

struct Error {
    std::string message;
};

// Parses, validates. All problems are reported together; the error message
// contains one bullet per issue with file:line where available.
std::expected<Config, Error> load(const std::filesystem::path& toml_file);

// Pure cross-field validation of an already-parsed Config.
std::expected<void, Error> validate(const Config&);

// Canonical TOML serialization (used by tools/cfg2toml and round-trip tests).
std::string to_toml(const Config&);

// Enum <-> string tables (exposed for tests and diagnostics).
std::string_view to_string(InviscidScheme);
std::string_view to_string(BcType);
std::string_view to_string(Edge);
std::optional<InviscidScheme> scheme_from_string(std::string_view);
std::optional<BcType> bc_from_string(std::string_view);
std::optional<Edge> edge_from_string(std::string_view);

}  // namespace ns::config
