#include "io/restart.hpp"

#include <cstring>
#include <fstream>

namespace ns::io {

namespace {
constexpr char kMagic[8] = {'N', 'S', 'R', 'S', 'T', '0', '0', '1'};
constexpr std::int32_t kVersion = 1;

std::size_t canvas_size(int ni, int nj, int ng) {
    return static_cast<std::size_t>(ni + 2 * ng) * static_cast<std::size_t>(nj + 2 * ng);
}

template <typename T>
bool write_pod(std::ostream& os, const T& v) {
    static_assert(std::is_trivially_copyable_v<T>);
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
    return os.good();
}

template <typename T>
bool read_pod(std::istream& is, T& v) {
    static_assert(std::is_trivially_copyable_v<T>);
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    return is.good();
}
}  // namespace

const RestartField* RestartData::find(std::string_view name) const {
    for (const auto& f : fields)
        if (f.name == name) return &f;
    return nullptr;
}

std::expected<void, Error> write_restart(const std::filesystem::path& path, int nx, int ny,
                                         int ng, FieldSpan fields) {
    std::ofstream os(path, std::ios::binary);
    if (!os) return std::unexpected(Error{std::format("cannot open '{}' for writing",
                                                      path.string())});

    os.write(kMagic, sizeof(kMagic));
    bool ok = write_pod(os, kVersion) && write_pod(os, nx) && write_pod(os, ny) &&
             write_pod(os, ng) && write_pod(os, static_cast<std::int32_t>(fields.size()));
    if (!ok) return std::unexpected(Error{"failed writing restart header"});

    for (const auto& [name, field] : fields) {
        const std::int32_t len = static_cast<std::int32_t>(name.size());
        if (!write_pod(os, len)) return std::unexpected(Error{"failed writing field name"});
        os.write(name.data(), name.size());
        if (!os) return std::unexpected(Error{"failed writing restart payload"});
        os.write(reinterpret_cast<const char*>(field->data()),
                 static_cast<std::streamsize>(field->size() * sizeof(double)));
        if (!os) return std::unexpected(Error{std::format("failed writing field '{}'", name)});
    }
    return {};
}

std::expected<RestartData, Error> read_restart(const std::filesystem::path& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) return std::unexpected(Error{std::format("cannot open '{}' for reading",
                                                      path.string())});

    char magic[sizeof(kMagic)];
    is.read(magic, sizeof(magic));
    if (!is || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0)
        return std::unexpected(Error{std::format("'{}' is not a solver restart file",
                                                 path.string())});

    RestartData data;
    std::int32_t version = 0, nf = 0;
    if (!read_pod(is, version) || version != kVersion)
        return std::unexpected(Error{"unsupported restart version"});
    if (!read_pod(is, data.nx) || !read_pod(is, data.ny) || !read_pod(is, data.ng) ||
        !read_pod(is, nf))
        return std::unexpected(Error{"truncated restart header"});

    for (std::int32_t i = 0; i < nf; ++i) {
        std::int32_t len = 0;
        if (!read_pod(is, len) || len < 0 || len > 4096)
            return std::unexpected(Error{"corrupt field name length"});
        RestartField f;
        f.name.resize(static_cast<std::size_t>(len));
        is.read(f.name.data(), len);
        if (!is) return std::unexpected(Error{"truncated field name"});

        f.values.resize(canvas_size(data.nx, data.ny, data.ng));
        is.read(reinterpret_cast<char*>(f.values.data()),
                static_cast<std::streamsize>(f.values.size() * sizeof(double)));
        if (!is) return std::unexpected(Error{std::format("truncated field '{}'", f.name)});
        data.fields.push_back(std::move(f));
    }
    return data;
}

}  // namespace ns::io
