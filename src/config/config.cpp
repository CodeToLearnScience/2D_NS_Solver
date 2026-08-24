#include "config/config.hpp"

#include <cmath>
#include <format>
#include <fstream>
#include <iomanip>
#include <set>
#include <span>
#include <sstream>
#include <utility>

#include <toml++/toml.hpp>

namespace ns::config {
namespace {

using toml::node;
using toml::table;

Error join_errors(std::vector<Error> errs, std::string_view header) {
    std::ostringstream out;
    out << header << "\n";
    for (const auto& e : errs) out << "  - " << e.message << "\n";
    return Error{out.str()};
}

// Numeric fields accept both TOML integers and floats (scaling = 1 must load).
std::optional<double> as_double(const node& n) {
    if (auto v = n.value<double>()) return v;
    if (auto v = n.value<int64_t>()) return static_cast<double>(*v);
    return std::nullopt;
}

template <typename T>
concept TomlScalar =
    std::same_as<T, std::string> || std::same_as<T, bool> ||
    std::same_as<T, double> || std::same_as<T, int>;

// Reads keys from one TOML table, accumulating errors with key paths, and
// remembering consumed keys so that unknown ones can be reported afterwards.
class TableReader {
public:
    explicit TableReader(const table* t, std::string path)
        : tbl_(t), path_(std::move(path)) {}

    template <TomlScalar T>
    std::optional<T> required(std::string_view key) {
        touch(key);
        const node* n = find(key);
        if (n == nullptr) {
            fail(std::format("missing required key '{}.{}'", path_, key));
            return std::nullopt;
        }
        return coerce<T>(*n, full(key));
    }

    template <TomlScalar T>
    std::optional<T> optional(std::string_view key) {
        touch(key);
        const node* n = find(key);
        if (n == nullptr) return std::nullopt;
        return coerce<T>(*n, full(key));
    }

    void touch(std::string_view key) { consumed_.emplace(key); }

    template <typename E>
    void required_enum(std::string_view key, std::string_view what,
                       const std::vector<std::pair<std::string_view, E>>& names, E& dest) {
        std::optional<E> got;
        enum_lookup(key, what, names, got);
        if (got) dest = *got;
    }

    template <typename E>
    void optional_enum(std::string_view key, std::string_view what,
                       const std::vector<std::pair<std::string_view, E>>& names,
                       std::optional<E>& dest) {
        enum_lookup(key, what, names, dest);
    }

    void fail(std::string msg) const { errors_.push_back(Error{std::move(msg)}); }

    void check_unknown_keys() const {
        if (tbl_ == nullptr) return;
        for (const auto& [k, v] : *tbl_) {
            if (!consumed_.contains(std::string{k}))
                fail(std::format("unknown key '{}.{}'", path_, std::string{k}));
        }
    }

    [[nodiscard]] bool ok() const { return errors_.empty(); }
    [[nodiscard]] std::vector<Error> take_errors() { return std::exchange(errors_, {}); }

    [[nodiscard]] std::string child_path(std::string_view key) const {
        return std::format("{}.{}", path_, key);
    }

private:
    const node* find(std::string_view key) const {
        if (tbl_ == nullptr) return nullptr;
        const node* n = tbl_->get(key);
        return (n != nullptr && !n->is_table()) ? n : nullptr;
    }

    template <TomlScalar T>
    std::optional<T> coerce(const node& n, const std::string& what) {
        if constexpr (std::same_as<T, double>) {
            if (auto v = as_double(n)) return v;
        } else if constexpr (std::same_as<T, int>) {
            if (auto v = n.value<int64_t>()) {
                return static_cast<int>(*v);
            }
        } else {
            if (auto v = n.value<T>()) return *v;
        }
        fail(std::format("{}: expected {}", what, expected_type_name<T>()));
        return std::nullopt;
    }

    template <TomlScalar T>
    static std::string_view expected_type_name() {
        if constexpr (std::same_as<T, std::string>) return "a string";
        else if constexpr (std::same_as<T, bool>) return "a boolean";
        else if constexpr (std::same_as<T, int>) return "an integer";
        else return "a number";
    }

    template <typename E>
    void enum_lookup(std::string_view key, std::string_view what,
                     const std::vector<std::pair<std::string_view, E>>& names,
                     std::optional<E>& dest) {
        auto s = optional<std::string>(key);
        if (!s) return;
        for (const auto& [name, value] : names) {
            if (*s == name) {
                dest = value;
                return;
            }
        }
        std::string allowed;
        for (const auto& [name, _] : names) {
            if (!allowed.empty()) allowed += ", ";
            allowed += std::format("'{}'", name);
        }
        fail(std::format("{}: invalid {} '{}' (expected one of {})", full(key), what, *s,
                         allowed));
    }

    std::string full(std::string_view key) const { return std::format("{}.{}", path_, key); }

    const table* tbl_;
    std::string path_;
    std::set<std::string> consumed_;
    mutable std::vector<Error> errors_;
};

inline const std::vector<std::pair<std::string_view, Equations>> kEquationsNames{
    {"euler", Equations::Euler},
    {"navier-stokes", Equations::NavierStokes}};
inline const std::vector<std::pair<std::string_view, Formulation>> kFormulationNames{
    {"nondimensional", Formulation::Nondimensional},
    {"dimensional", Formulation::Dimensional}};
inline const std::vector<std::pair<std::string_view, FlowType>> kFlowTypeNames{
    {"steady", FlowType::Steady}, {"unsteady", FlowType::Unsteady}};
inline const std::vector<std::pair<std::string_view, TimeMethod>> kTimeMethodNames{
    {"forward-euler", TimeMethod::ForwardEuler},
    {"ssprk2", TimeMethod::SSPRK2},
    {"ssprk3", TimeMethod::SSPRK3}};
inline const std::vector<std::pair<std::string_view, SpatialOrder>> kOrderNames{
    {"first", SpatialOrder::First}, {"second", SpatialOrder::Second}};
inline const std::vector<std::pair<std::string_view, ViscousMethod>> kViscousNames{
    {"green-gauss", ViscousMethod::GreenGauss},
    {"least-squares", ViscousMethod::LeastSquares}};

template <typename E>
std::vector<std::pair<std::string_view, E>> full_enum_table(
    std::string_view (*name_of)(E), int last) {
    std::vector<std::pair<std::string_view, E>> t;
    for (int i = 0; i <= last; ++i) {
        auto v = static_cast<E>(i);
        t.emplace_back(name_of(v), v);
    }
    return t;
}

std::string_view equations_name(Equations e) {
    return e == Equations::Euler ? "euler" : "navier-stokes";
}

void put_d(std::ostream& os, double v) {
    // Shortest representation that still round-trips exactly.
    std::string s = std::isfinite(v) ? std::format("{}", v) : "nan";
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos &&
        s.find('n') == std::string::npos && s.find('i') == std::string::npos) {
        s += ".0";  // keep it a TOML float so it loads as a number
    }
    os << s;
}

std::expected<Config, Error> parse(const table& root) {
    Config cfg;

    {  // [case]
        TableReader r{root.get_as<table>("case"), "case"};
        if (auto v = r.required<std::string>("name")) cfg.case_name = std::move(*v);
        r.required_enum("equations", "equations type", kEquationsNames, cfg.equations);
        r.required_enum("formulation", "formulation", kFormulationNames, cfg.formulation);
        r.required_enum("flow", "flow type", kFlowTypeNames, cfg.flow);
        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[case]"));
    }

    {  // [grid]
        TableReader r{root.get_as<table>("grid"), "grid"};
        if (auto v = r.required<std::string>("file")) cfg.grid.file = std::move(*v);
        if (auto v = r.required<int>("nx")) cfg.grid.nx = *v;
        if (auto v = r.required<int>("ny")) cfg.grid.ny = *v;
        if (auto v = r.optional<double>("scaling")) cfg.grid.scaling = *v;
        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[grid]"));
    }

    {  // [numerics]
        TableReader r{root.get_as<table>("numerics"), "numerics"};
        r.required_enum("inviscid_scheme", "inviscid scheme",
                        full_enum_table<InviscidScheme>(to_string,
                                                        static_cast<int>(InviscidScheme::NKFDS_MOVERS_2O)),
                        cfg.numerics.inviscid_scheme);
        r.required_enum("order", "spatial order", kOrderNames, cfg.numerics.order);
        r.required_enum("time_method", "time-stepping method", kTimeMethodNames,
                        cfg.numerics.time_method);
        if (auto v = r.required<double>("cfl")) cfg.numerics.cfl = *v;
        std::optional<ViscousMethod> vm;
        r.optional_enum("viscous_method", "viscous discretization", kViscousNames, vm);
        if (vm) cfg.numerics.viscous_method = *vm;
        if (auto v = r.optional<double>("residual_smoothing_eps"))
            cfg.numerics.residual_smoothing_eps = *v;
        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[numerics]"));
    }

    {  // [physics] + [physics.reference_state]
        TableReader r{root.get_as<table>("physics"), "physics"};
        if (auto v = r.optional<double>("gamma")) cfg.physics.gamma = *v;
        if (auto v = r.optional<double>("prandtl")) cfg.physics.prandtl = *v;
        if (auto v = r.optional<double>("gas_constant")) cfg.physics.gas_constant = *v;
        if (auto v = r.optional<double>("Re_inf")) cfg.physics.re_inf = *v;
        if (auto v = r.optional<double>("Mach_inf")) cfg.physics.mach_inf = *v;
        if (auto v = r.optional<double>("alpha_deg")) cfg.physics.alpha_deg = *v;
        if (auto v = r.optional<double>("p_inf")) cfg.physics.p_inf = *v;
        if (auto v = r.optional<double>("T_inf")) cfg.physics.t_inf = *v;
        if (auto v = r.optional<double>("reference_length"))
            cfg.physics.reference_length = *v;

        PhysicsConfig::ReferenceState ref{};
        bool have_ref = false;
        if (auto rs = root.at_path("physics.reference_state"); rs && rs.is_table()) {
            r.touch("reference_state");
            TableReader rr{rs.as_table(), r.child_path("reference_state")};
            if (auto v = rr.optional<double>("pressure")) { ref.pressure = *v; have_ref = true; }
            if (auto v = rr.optional<double>("temperature")) { ref.temperature = *v; have_ref = true; }
            if (auto v = rr.optional<double>("density")) { ref.density = *v; have_ref = true; }
            if (auto v = rr.optional<double>("velocity")) { ref.velocity = *v; have_ref = true; }
            rr.check_unknown_keys();
            if (!rr.ok())
                return std::unexpected(
                    join_errors(rr.take_errors(), "[physics.reference_state]"));
        } else if (rs && rs.is_value()) {
            r.fail("[physics.reference_state]: expected a table");
        }
        if (have_ref) cfg.physics.reference_state = ref;

        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[physics]"));
    }

    {  // [run]
        TableReader r{root.get_as<table>("run"), "run"};
        if (auto v = r.required<int>("max_iterations")) cfg.run.max_iterations = *v;
        if (auto v = r.optional<double>("total_time")) cfg.run.total_time = *v;
        if (auto v = r.optional<bool>("restart")) cfg.run.restart = *v;
        if (auto v = r.optional<std::string>("restart_file"))
            cfg.run.restart_file = std::move(*v);
        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[run]"));
    }

    {  // [output]
        TableReader r{root.get_as<table>("output"), "output"};
        if (auto v = r.optional<std::string>("directory"))
            cfg.output.directory = std::move(*v);
        if (auto v = r.optional<int>("display_frequency"))
            cfg.output.display_frequency = *v;
        if (auto v = r.optional<int>("write_frequency"))
            cfg.output.write_frequency = *v;
        if (auto fmt_node = root.at_path("output.formats")) {
            r.touch("formats");
            if (const toml::array* arr = fmt_node.as_array()) {
                cfg.output.formats.clear();
                for (const auto& e : *arr) {
                    if (auto s = e.value<std::string>())
                        cfg.output.formats.push_back(*s);
                    else {
                        r.fail("output.formats: entries must be strings");
                        break;
                    }
                }
            } else {
                r.fail("output.formats: expected an array of strings");
            }
        }
        r.check_unknown_keys();
        if (!r.ok()) return std::unexpected(join_errors(r.take_errors(), "[output]"));
    }

    {  // [[boundary]]
        if (auto segs = root.at_path("boundary"); segs && segs.is_array()) {
            int idx = -1;
            for (const auto& elem : *segs.as_array()) {
                ++idx;
                const table* st = elem.as_table();
                if (st == nullptr)
                    return std::unexpected(
                        Error{std::format("[[boundary]] #{}: expected a table", idx)});
                BoundarySegment seg{};
                TableReader r{st, std::format("boundary[{}]", idx)};
                r.required_enum("edge", "edge name",
                                full_enum_table<Edge>(to_string, static_cast<int>(Edge::YMax)),
                                seg.edge);
                if (auto v = r.required<int>("start")) seg.start = *v;
                if (auto v = r.required<int>("end")) seg.end = *v;
                r.required_enum("type", "boundary type",
                                full_enum_table<BcType>(to_string, static_cast<int>(BcType::Cut)),
                                seg.type);

                if (auto state = st->at_path("state"); state && state.is_array()) {
                    r.touch("state");
                    std::array<double, 4> vals{};
                    std::size_t c = 0;
                    bool all_num = true;
                    for (const auto& e : *state.as_array()) {
                        if (c >= vals.size()) { all_num = false; break; }
                        if (auto d = as_double(e)) vals[c++] = *d;
                        else { all_num = false; break; }
                    }
                    if (all_num && c == 4) {
                        seg.prescribed_state = vals;
                    } else {
                        r.fail(std::format("boundary[{}].state: expected 4 numbers "
                                           "(rho, u, v, p)",
                                           idx));
                    }
                } else if (state && state.is_value()) {
                    r.touch("state");
                    r.fail(std::format("boundary[{}].state: expected an array of 4 numbers "
                                       "(rho, u, v, p)",
                                       idx));
                }

                if (auto v = r.optional<double>("wall_temperature"))
                    seg.wall_temperature = *v;

                r.check_unknown_keys();
                if (!r.ok())
                    return std::unexpected(join_errors(
                        r.take_errors(), std::format("[[boundary]] #{}", idx)));
                cfg.boundaries.push_back(std::move(seg));
            }
        }
    }

    return cfg;
}

}  // namespace

std::expected<void, Error> validate(const Config& c) {
    std::vector<Error> errs;
    auto add = [&errs](std::string msg) { errs.push_back(Error{std::move(msg)}); };

    if (c.case_name.empty()) add("case.name must not be empty");
    if (c.grid.file.empty()) add("grid.file must not be empty");
    if (c.grid.nx < 1) add("grid.nx must be >= 1");
    if (c.grid.ny < 1) add("grid.ny must be >= 1");
    if (!(c.grid.scaling > 0.0)) add("grid.scaling must be > 0");

    if (!(c.numerics.cfl > 0.0)) add("numerics.cfl must be > 0");
    if (c.numerics.cfl > 10.0) add("numerics.cfl > 10 is implausible; refusing to run");

    if (c.output.display_frequency < 1) add("output.display_frequency must be >= 1");
    if (c.output.write_frequency < 1) add("output.write_frequency must be >= 1");
    if (c.output.directory.empty()) add("output.directory must not be empty");

    if (c.run.max_iterations < 1 && c.flow == FlowType::Steady)
        add("run.max_iterations must be >= 1 for steady runs");
    if (c.run.total_time <= 0.0 && c.flow == FlowType::Unsteady)
        add("run.total_time must be > 0 for unsteady runs");
    if (c.run.restart && c.run.restart_file.empty())
        add("run.restart_file is required when run.restart = true");

    if (!(c.physics.gamma > 1.0)) add("physics.gamma must be > 1");
    if (!(c.physics.prandtl > 0.0)) add("physics.prandtl must be > 0");
    if (!(c.physics.gas_constant > 0.0)) add("physics.gas_constant must be > 0");
    if (!(c.physics.mach_inf > 0.0)) add("physics.Mach_inf must be > 0");
    if (!(c.physics.re_inf > 0.0)) add("physics.Re_inf must be > 0");
    if (c.formulation == Formulation::Dimensional) {
        // Only consumed by the dimensional initialization path.
        if (!(c.physics.p_inf > 0.0)) add("physics.p_inf must be > 0 for dimensional runs");
        if (!(c.physics.t_inf > 0.0)) add("physics.T_inf must be > 0 for dimensional runs");
    }
    if (!(c.physics.reference_length > 0.0))
        add("physics.reference_length must be > 0");

    if (c.formulation == Formulation::Nondimensional && !c.physics.reference_state.has_value())
        add("physics.reference_state is required for nondimensional runs");
    if (c.physics.reference_state) {
        const auto& r = *c.physics.reference_state;
        if (!(r.pressure > 0.0) || !(r.temperature > 0.0) || !(r.density > 0.0) ||
            !(r.velocity > 0.0))
            add("physics.reference_state values must all be > 0");
    }

    const std::set<std::string> valid_formats{"solution", "vtk", "surface", "residual"};
    if (c.output.formats.empty()) add("output.formats must not be empty");
    for (const auto& f : c.output.formats) {
        if (!valid_formats.contains(f))
            add(std::format("output.formats: unknown format '{}' (valid: solution, vtk, "
                            "surface, residual)",
                            f));
    }

    if (c.boundaries.empty()) add("at least one [[boundary]] segment is required");

    for (std::size_t i = 0; i < c.boundaries.size(); ++i) {
        const auto& b = c.boundaries[i];
        const std::string tag = std::format("boundary[{}]", i);
        if (b.start < 2) add(std::format("{}: start must be >= 2 (legacy node indexing)", tag));
        if (b.end < b.start) add(std::format("{}: end must be >= start", tag));

        const int limit =
            (b.edge == Edge::XMin || b.edge == Edge::XMax) ? c.grid.ny + 1 : c.grid.nx + 1;
        if (b.end > limit)
            add(std::format("{}: end={} exceeds the last node ({}) for this edge", tag, b.end,
                            limit));

        switch (b.type) {
            case BcType::PrescribedInflow:
                // Optional 'state': when absent the free-stream state is used,
                // matching legacy pres_input_flag = 0 semantics.
                break;
            case BcType::NoSlipIsothermalWall:
                if (!b.wall_temperature.has_value() || !(*b.wall_temperature > 0.0))
                    add(std::format(
                        "{}: no-slip-isothermal-wall requires wall_temperature > 0", tag));
                break;
            default:
                if (b.prescribed_state.has_value())
                    add(std::format("{}: 'state' is only valid for prescribed-inflow", tag));
                if (b.wall_temperature.has_value())
                    add(std::format(
                        "{}: wall_temperature is only valid for no-slip-isothermal-wall",
                        tag));
                break;
        }
    }

    if (errs.empty()) return {};
    return std::unexpected(
        join_errors(errs, std::format("invalid configuration '{}':", c.case_name)));
}

std::expected<Config, Error> load(const std::filesystem::path& toml_file) {
    std::error_code ec;
    if (!std::filesystem::is_regular_file(toml_file, ec))
        return std::unexpected(
            Error{std::format("config file not found: '{}'", toml_file.string())});

    toml::table root;
    try {
        root = toml::parse_file(toml_file.string());
    } catch (const toml::parse_error& e) {
        std::ostringstream out;
        out << "TOML parse error in '" << toml_file.string() << "':\n"
            << e.description() << "\n";
        return std::unexpected(Error{out.str()});
    }

    auto parsed = parse(root);
    if (!parsed) return std::unexpected(parsed.error());

    if (auto v = validate(parsed.value()); !v) return std::unexpected(v.error());
    return std::move(parsed.value());
}

// ---------------------------------------------------------------------------
// enum <-> string
// ---------------------------------------------------------------------------
std::string_view to_string(InviscidScheme s) {
    switch (s) {
        case InviscidScheme::LLF: return "llf";
        case InviscidScheme::MOVERS: return "movers";
        case InviscidScheme::MOVERS_H1: return "movers-h1";
        case InviscidScheme::MOVERS_LE1: return "movers-le1";
        case InviscidScheme::ROE_TV: return "roe-tv";
        case InviscidScheme::KFDS: return "kfds";
        case InviscidScheme::NKFDS_MOVERS: return "nkfds-movers";
        case InviscidScheme::LLF_2O: return "llf-2o";
        case InviscidScheme::MOVERS_H2: return "movers-h2";
        case InviscidScheme::MOVERS_LE2: return "movers-le2";
        case InviscidScheme::ROE_2O: return "roe-2o";
        case InviscidScheme::ROE_TV_2O: return "roe-tv-2o";
        case InviscidScheme::KFDS_2O: return "kfds-2o";
        case InviscidScheme::ECCS: return "eccs";
        case InviscidScheme::MOVERS_NWSC: return "movers-nwsc";
        case InviscidScheme::NKFDS_MOVERS_2O: return "nkfds-movers-2o";
    }
    return "?";
}

std::string_view to_string(BcType b) {
    switch (b) {
        case BcType::PrescribedInflow: return "prescribed-inflow";
        case BcType::Transmissive: return "transmissive";
        case BcType::SlipWall: return "slip-wall";
        case BcType::NoSlipAdiabaticWall: return "no-slip-adiabatic-wall";
        case BcType::NoSlipIsothermalWall: return "no-slip-isothermal-wall";
        case BcType::Farfield: return "farfield";
        case BcType::Symmetry: return "symmetry";
        case BcType::Cut: return "cut";
    }
    return "?";
}

std::string_view to_string(Edge e) {
    switch (e) {
        case Edge::XMin: return "imin";
        case Edge::XMax: return "imax";
        case Edge::YMin: return "jmin";
        case Edge::YMax: return "jmax";
    }
    return "?";
}

std::optional<InviscidScheme> scheme_from_string(std::string_view s) {
    for (int i = 0; i <= static_cast<int>(InviscidScheme::NKFDS_MOVERS_2O); ++i) {
        auto v = static_cast<InviscidScheme>(i);
        if (to_string(v) == s) return v;
    }
    return std::nullopt;
}

std::optional<BcType> bc_from_string(std::string_view s) {
    for (int i = 0; i <= static_cast<int>(BcType::Cut); ++i) {
        auto v = static_cast<BcType>(i);
        if (to_string(v) == s) return v;
    }
    return std::nullopt;
}

std::optional<Edge> edge_from_string(std::string_view s) {
    for (int i = 0; i <= static_cast<int>(Edge::YMax); ++i) {
        auto v = static_cast<Edge>(i);
        if (to_string(v) == s) return v;
    }
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// canonical TOML emission
// ---------------------------------------------------------------------------
std::string to_toml(const Config& c) {
    std::ostringstream os;
    os << "# Generated by ns::config::to_toml -- canonical schema\n\n";

    os << "[case]\n";
    os << "name = \"" << c.case_name << "\"\n";
    os << "equations = \"" << equations_name(c.equations) << "\"\n";
    os << "formulation = \""
       << (c.formulation == Formulation::Nondimensional ? "nondimensional" : "dimensional")
       << "\"\n";
    os << "flow = \"" << (c.flow == FlowType::Steady ? "steady" : "unsteady") << "\"\n\n";

    os << "[grid]\n";
    os << "file = \"" << c.grid.file << "\"\n";
    os << "nx = " << c.grid.nx << "\n";
    os << "ny = " << c.grid.ny << "\n";
    os << "scaling = "; put_d(os, c.grid.scaling); os << "\n\n";

    os << "[numerics]\n";
    os << "inviscid_scheme = \"" << to_string(c.numerics.inviscid_scheme) << "\"\n";
    os << "order = \"" << (c.numerics.order == SpatialOrder::First ? "first" : "second")
       << "\"\n";
    std::string_view tm = "forward-euler";
    switch (c.numerics.time_method) {
        case TimeMethod::SSPRK2: tm = "ssprk2"; break;
        case TimeMethod::SSPRK3: tm = "ssprk3"; break;
        default: break;
    }
    os << "time_method = \"" << tm << "\"\n";
    os << "cfl = "; put_d(os, c.numerics.cfl); os << "\n";
    os << "viscous_method = \""
       << (c.numerics.viscous_method == ViscousMethod::GreenGauss ? "green-gauss"
                                                                  : "least-squares")
       << "\"\n";
    if (c.numerics.residual_smoothing_eps != 0.0) {
        os << "residual_smoothing_eps = ";
        put_d(os, c.numerics.residual_smoothing_eps);
        os << "\n";
    }
    os << "\n";

    os << "[physics]\n";
    if (c.physics.gamma != 1.4) { os << "gamma = "; put_d(os, c.physics.gamma); os << "\n"; }
    if (c.physics.prandtl != 0.72) { os << "prandtl = "; put_d(os, c.physics.prandtl); os << "\n"; }
    if (c.physics.gas_constant != 286.92) {
        os << "gas_constant = "; put_d(os, c.physics.gas_constant); os << "\n";
    }
    os << "Re_inf = "; put_d(os, c.physics.re_inf); os << "\n";
    os << "Mach_inf = "; put_d(os, c.physics.mach_inf); os << "\n";
    os << "alpha_deg = "; put_d(os, c.physics.alpha_deg); os << "\n";
    os << "p_inf = "; put_d(os, c.physics.p_inf); os << "\n";
    os << "T_inf = "; put_d(os, c.physics.t_inf); os << "\n";
    if (c.physics.reference_length != 1.0) {
        os << "reference_length = ";
        put_d(os, c.physics.reference_length);
        os << "\n";
    }
    if (c.physics.reference_state) {
        os << "\n[physics.reference_state]\n";
        os << "pressure = ";    put_d(os, c.physics.reference_state->pressure);    os << "\n";
        os << "temperature = "; put_d(os, c.physics.reference_state->temperature); os << "\n";
        os << "density = ";     put_d(os, c.physics.reference_state->density);     os << "\n";
        os << "velocity = ";    put_d(os, c.physics.reference_state->velocity);    os << "\n";
    }
    os << "\n";

    os << "[run]\n";
    os << "max_iterations = " << c.run.max_iterations << "\n";
    if (c.run.total_time > 0.0) {
        os << "total_time = ";
        put_d(os, c.run.total_time);
        os << "\n";
    }
    if (c.run.restart) {
        os << "restart = true\n";
        os << "restart_file = \"" << c.run.restart_file << "\"\n";
    }
    os << "\n";

    os << "[output]\n";
    os << "directory = \"" << c.output.directory << "\"\n";
    os << "display_frequency = " << c.output.display_frequency << "\n";
    os << "write_frequency = " << c.output.write_frequency << "\n";
    os << "formats = [";
    for (std::size_t i = 0; i < c.output.formats.size(); ++i) {
        if (i != 0) os << ", ";
        os << '"' << c.output.formats[i] << '"';
    }
    os << "]\n";

    if (!c.boundaries.empty()) {
        os << "\n";
        for (const auto& b : c.boundaries) {
            os << "[[boundary]]\n";
            os << "edge = \"" << to_string(b.edge) << "\"\n";
            os << "start = " << b.start << "\n";
            os << "end = " << b.end << "\n";
            os << "type = \"" << to_string(b.type) << "\"\n";
            if (b.prescribed_state) {
                os << "state = [";
                for (int k = 0; k < 4; ++k) {
                    if (k != 0) os << ", ";
                    put_d(os, (*b.prescribed_state)[static_cast<std::size_t>(k)]);
                }
                os << "]\n";
            }
            if (b.wall_temperature) {
                os << "wall_temperature = ";
                put_d(os, *b.wall_temperature);
                os << "\n";
            }
            os << "\n";
        }
    }

    return os.str();
}

}  // namespace ns::config
