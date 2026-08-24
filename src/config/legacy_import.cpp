#include "config/legacy_import.hpp"

#include <array>
#include <charconv>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace ns::config {
namespace {

// Whitespace tokenizer over a legacy config file. Comment text ('#' to end of
// line) is stripped before tokenizing, so header/comment lines vanish and the
// remaining tokens are exactly the positional values the old reader consumed.
class TokenStream {
public:
    explicit TokenStream(const std::filesystem::path& p) {
        std::ifstream in(p);
        if (!in) {
            failed_ = true;
            return;
        }
        std::ostringstream raw;
        raw << in.rdbuf();
        std::istringstream line_stream(raw.str());
        std::string line;
        while (std::getline(line_stream, line)) {
            if (auto pos = line.find('#'); pos != std::string::npos) line.resize(pos);
            buffer_ += line;
            buffer_ += ' ';
        }
        tokens_ = split(buffer_);
    }

    [[nodiscard]] bool ok() const { return !failed_; }
    [[nodiscard]] std::size_t size() const { return tokens_.size(); }
    [[nodiscard]] std::size_t pos() const { return pos_; }

    std::optional<std::string> next_string() {
        if (pos_ >= tokens_.size()) return std::nullopt;
        return tokens_[pos_++];
    }

    template <typename T>
    std::optional<T> next_number() {
        auto s = next_string();
        if (!s) return std::nullopt;
        T v{};
        const char* begin = s->data();
        const char* end = s->data() + s->size();
        if constexpr (std::is_integral_v<T>) {
            auto [ptr, ec] = std::from_chars(begin, end, v);
            if (ec != std::errc{} || ptr != end) return std::nullopt;
        } else {
            // from_chars for floating point is fine on GCC; strtod fallback not needed.
            auto [ptr, ec] = std::from_chars(begin, end, v);
            if (ec != std::errc{} || ptr != end) return std::nullopt;
        }
        return v;
    }

private:
    static std::vector<std::string> split(const std::string& s) {
        std::vector<std::string> out;
        std::istringstream ss(s);
        std::string tok;
        while (ss >> tok) out.push_back(tok);
        return out;
    }

    bool failed_ = false;
    std::string buffer_;
    std::vector<std::string> tokens_;
    std::size_t pos_ = 0;
};

struct TopoSegment {
    int pres_input_flag = 0;
    BcType type{};
    Edge edge{};
    int start = 0;
    int end = 0;
    std::optional<std::array<double, 4>> state;
    std::optional<double> wall_temperature;
};

Error fail_ctx(const std::filesystem::path& p, std::string_view what) {
    return Error{std::format("legacy import of '{}': {}", p.string(), what)};
}

std::expected<Edge, Error> edge_from_legacy_ind(int ind, const std::filesystem::path& p) {
    switch (ind) {  // 1=bottom, 2=right, 3=top, 4=left (bc_transmitive.cpp diagnostics)
        case 1: return Edge::YMin;
        case 2: return Edge::XMax;
        case 3: return Edge::YMax;
        case 4: return Edge::XMin;
        default:
            return std::unexpected(
                fail_ctx(p, std::format("invalid boundary index {} (expected 1..4)", ind)));
    }
}

std::expected<std::vector<TopoSegment>, Error> read_topology(
    const std::filesystem::path& topo_path) {
    TokenStream ts(topo_path);
    if (!ts.ok()) return std::unexpected(fail_ctx(topo_path, "file could not be opened"));

    std::vector<TopoSegment> segs;
    try {
        const int n_seg = ts.next_number<int>().value_or(-1);
        const int n_presc = ts.next_number<int>().value_or(0);
        if (n_seg < 0) return std::unexpected(fail_ctx(topo_path, "missing segment count"));

        int presc_seen = 0;
        segs.reserve(static_cast<std::size_t>(n_seg));
        for (int i = 0; i < n_seg; ++i) {
            TopoSegment s{};
            s.pres_input_flag = ts.next_number<int>().value_or(-1);
            const int bc_flag = ts.next_number<int>().value_or(-1);
            const int bound_ind = ts.next_number<int>().value_or(-1);

            // bound_cell (ghost-row index) is ignored: derivable from edge.
            if (!ts.next_number<int>().has_value())
                return std::unexpected(fail_ctx(topo_path, "truncated segment table"));

            s.start = ts.next_number<int>().value_or(-1);
            s.end = ts.next_number<int>().value_or(-1);

            switch (bc_flag) {
                case 10: s.type = BcType::PrescribedInflow; break;
                case 20: s.type = BcType::Transmissive; break;
                case 30: s.type = BcType::SlipWall; break;
                case 40: s.type = BcType::NoSlipAdiabaticWall; break;
                case 50: s.type = BcType::NoSlipIsothermalWall; break;
                case 60: s.type = BcType::Farfield; break;
                case 70: s.type = BcType::Symmetry; break;
                case 80: s.type = BcType::Cut; break;
                default:
                    return std::unexpected(fail_ctx(
                        topo_path,
                        std::format("segment {}: unknown bc flag {}", i, bc_flag)));
            }
            auto edge_res = edge_from_legacy_ind(bound_ind, topo_path);
            if (!edge_res) return std::unexpected(edge_res.error());
            s.edge = *edge_res;

            if (s.pres_input_flag == 1) {
                std::array<double, 4> st{};
                for (double& v : st)
                    if (auto x = ts.next_number<double>()) v = *x;
                    else return std::unexpected(fail_ctx(topo_path, "truncated prescribed state"));
                s.state = st;
                ++presc_seen;
            } else if (bc_flag == 50) {
                if (auto t = ts.next_number<double>()) s.wall_temperature = *t;
                else return std::unexpected(fail_ctx(topo_path, "truncated isothermal wall"));
            }
            segs.push_back(std::move(s));
        }
        if (presc_seen != n_presc)
            return std::unexpected(fail_ctx(
                topo_path,
                std::format("declared {} prescribed segments but found {}", n_presc,
                            presc_seen)));
    } catch (...) {
        return std::unexpected(fail_ctx(topo_path, "malformed topology file"));
    }
    return segs;
}

}  // namespace

std::expected<Config, Error> import_legacy(const std::filesystem::path& legacy_cfg) {
    std::error_code ec;
    if (!std::filesystem::is_regular_file(legacy_cfg, ec))
        return std::unexpected(Error{std::format("legacy config not found: '{}'",
                                                 legacy_cfg.string())});

    TokenStream ts(legacy_cfg);
    if (!ts.ok()) return std::unexpected(fail_ctx(legacy_cfg, "file could not be opened"));

    Config c;
    c.output.directory = "output";

    auto need_str = [&](std::string_view field) -> std::expected<std::string, Error> {
        auto v = ts.next_string();
        if (!v) return std::unexpected(fail_ctx(legacy_cfg, std::format("missing '{}'", field)));
        return std::move(*v);
    };
    auto need_d = [&](std::string_view field) -> std::expected<double, Error> {
        auto v = ts.next_number<double>();
        if (!v) return std::unexpected(fail_ctx(legacy_cfg, std::format("missing '{}'", field)));
        return *v;
    };
    auto need_i = [&](std::string_view field) -> std::expected<int, Error> {
        auto v = ts.next_number<int>();
        if (!v) return std::unexpected(fail_ctx(legacy_cfg, std::format("missing '{}'", field)));
        return *v;
    };

    // ---- positional block, order fixed by read_solver_input.cpp ------------
    if (auto v = need_str("test_case")) c.case_name = std::move(*v); else return std::unexpected(v.error());
    const int eqn_type = need_i("eqn_type").value_or(-1);
    const int dimen = need_i("dimen").value_or(-1);
    const int flow_type = need_i("flow_type").value_or(-1);
    const int time_accuracy = need_i("time_accuracy").value_or(-1);
    const int inviscid_scheme = need_i("inviscid_scheme").value_or(-1);
    const int space_accuracy = need_i("space_accuracy").value_or(-1);
    const int visc_method = need_i("visc_method").value_or(0);
    if (auto v = need_d("cfl")) c.numerics.cfl = *v; else return std::unexpected(v.error());
    double scaling = need_d("scaling_factor").value_or(1.0);
    const int restart = need_i("restart").value_or(0);
    std::string mesh_file = need_str("mesh_file").value_or("");
    const int nx = need_i("nx").value_or(0);
    const int ny = need_i("ny").value_or(0);
    std::string mesh_top_file = need_str("mesh_top_file").value_or("");
    std::string restart_file = need_str("restart_file").value_or("");
    if (auto v = need_i("disp_freq")) c.output.display_frequency = *v; else return std::unexpected(v.error());
    if (auto v = need_i("outp_freq")) c.output.write_frequency = *v; else return std::unexpected(v.error());
    const int max_iter = need_i("max_iter").value_or(0);
    const double tot_time = need_d("tot_time").value_or(0.0);

    c.physics.re_inf = need_d("Re_inf").value_or(0.0);
    c.physics.mach_inf = need_d("Mach_inf").value_or(0.0);
    c.physics.reference_length = need_d("Lref").value_or(1.0);
    c.physics.alpha_deg = need_d("alpha").value_or(0.0);
    c.physics.p_inf = need_d("p_inf").value_or(0.0);
    c.physics.t_inf = need_d("T_inf").value_or(0.0);

    PhysicsConfig::ReferenceState ref{};
    bool have_ref = false;
    if (dimen == 0) {
        if (auto v = need_d("pref")) { ref.pressure = *v; have_ref = true; } else return std::unexpected(v.error());
        if (auto v = need_d("Tref")) { ref.temperature = *v; have_ref = true; } else return std::unexpected(v.error());
        if (auto v = need_d("rhoref")) { ref.density = *v; have_ref = true; } else return std::unexpected(v.error());
        if (auto v = need_d("velref")) { ref.velocity = *v; have_ref = true; } else return std::unexpected(v.error());
    }

    // ---- translate enums ----------------------------------------------------
    switch (eqn_type) {
        case 0: c.equations = Equations::Euler; break;
        case 1: c.equations = Equations::NavierStokes; break;
        default: return std::unexpected(fail_ctx(legacy_cfg, std::format("bad eqn_type {}", eqn_type)));
    }
    switch (dimen) {
        case 0: c.formulation = Formulation::Nondimensional; break;
        case 1: c.formulation = Formulation::Dimensional; break;
        default: return std::unexpected(fail_ctx(legacy_cfg, std::format("bad dimen {}", dimen)));
    }
    switch (flow_type) {
        case 0: c.flow = FlowType::Steady; break;
        case 1: c.flow = FlowType::Unsteady; break;
        default: return std::unexpected(fail_ctx(legacy_cfg, std::format("bad flow_type {}", flow_type)));
    }
    switch (time_accuracy) {
        case 0: c.numerics.time_method = TimeMethod::ForwardEuler; break;
        case 1: c.numerics.time_method = TimeMethod::SSPRK2; break;
        case 2: c.numerics.time_method = TimeMethod::SSPRK3; break;
        default:
            return std::unexpected(
                fail_ctx(legacy_cfg, std::format("bad time_accuracy {}", time_accuracy)));
    }
    switch (space_accuracy) {
        case 0: c.numerics.order = SpatialOrder::First; break;
        case 1: c.numerics.order = SpatialOrder::Second; break;
        default: return std::unexpected(fail_ctx(legacy_cfg, std::format("bad space_accuracy {}", space_accuracy)));
    }
    c.run.max_iterations = max_iter;
    c.run.total_time = tot_time;
    if (restart == 1) {
        c.run.restart = true;
        c.run.restart_file = restart_file;
    }
    c.grid.file = mesh_file;
    c.grid.nx = nx;
    c.grid.ny = ny;
    c.grid.scaling = scaling;

    // inviscid scheme codes -> enumerators (see solver.cpp dispatch)
    using S = InviscidScheme;
    struct BadScheme {};
    try {
        if (space_accuracy == 0) {
            switch (inviscid_scheme) {
                case 0: c.numerics.inviscid_scheme = S::LLF; break;
                case 1: c.numerics.inviscid_scheme = S::MOVERS; break;
                case 2: c.numerics.inviscid_scheme = S::MOVERS_H1; break;
                case 3: c.numerics.inviscid_scheme = S::MOVERS_LE1; break;
                case 5: c.numerics.inviscid_scheme = S::ROE_TV; break;
                case 6: c.numerics.inviscid_scheme = S::KFDS; break;
                case 7: c.numerics.inviscid_scheme = S::NKFDS_MOVERS; break;
                default: throw BadScheme{};
            }
        } else {
            switch (inviscid_scheme) {
                case 20: c.numerics.inviscid_scheme = S::LLF_2O; break;
                case 22: c.numerics.inviscid_scheme = S::MOVERS_H2; break;
                case 23: c.numerics.inviscid_scheme = S::MOVERS_LE2; break;
                case 24: c.numerics.inviscid_scheme = S::ROE_2O; break;
                case 25: c.numerics.inviscid_scheme = S::ROE_TV_2O; break;
                case 26: c.numerics.inviscid_scheme = S::KFDS_2O; break;
                case 27: c.numerics.inviscid_scheme = S::ECCS; break;
                case 28: c.numerics.inviscid_scheme = S::MOVERS; break;  // sic: legacy calls 1st-order MOVERS
                case 29: c.numerics.inviscid_scheme = S::MOVERS_NWSC; break;
                case 30: c.numerics.inviscid_scheme = S::NKFDS_MOVERS_2O; break;
                default: throw BadScheme{};
            }
        }
    } catch (const BadScheme&) {
        return std::unexpected(fail_ctx(
            legacy_cfg,
            std::format("unsupported inviscid_scheme {} at space_accuracy {}", inviscid_scheme,
                        space_accuracy)));
    }

    switch (visc_method) {
        case 0: c.numerics.viscous_method = ViscousMethod::GreenGauss; break;
        case 1: c.numerics.viscous_method = ViscousMethod::LeastSquares; break;
        default: return std::unexpected(fail_ctx(legacy_cfg, std::format("bad visc_method {}", visc_method)));
    }

    if (have_ref) c.physics.reference_state = ref;

    // ---- topology ------------------------------------------------------------
    const std::filesystem::path cfg_dir =
        std::filesystem::absolute(legacy_cfg, ec).parent_path();
    auto segs = read_topology(cfg_dir / mesh_top_file);
    if (!segs) return std::unexpected(segs.error());

    for (const auto& s : *segs) {
        BoundarySegment b;
        b.edge = s.edge;
        b.start = s.start;
        b.end = s.end;
        b.type = s.type;
        b.prescribed_state = s.state;
        b.wall_temperature = s.wall_temperature;
        c.boundaries.push_back(b);
    }

    // Legacy files never carry these; defaults already set on Config.

    if (auto v = validate(c); !v) return std::unexpected(v.error());
    return c;
}

}  // namespace ns::config
