#pragma once

// Binary restart dump/load for field data.
//
// Format (host endianness, fixed-width types):
//   char[8]   magic "NSRST001"
//   int32     version = 1
//   int32     nx, ny, ng
//   int32     number of fields
//   per field: int32 name length | name bytes | ni*nj_total doubles (raw canvas)

#include <expected>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "core/error.hpp"
#include "fields/field.hpp"

namespace ns::io {

struct RestartField {
    std::string name;
    std::vector<double> values;  // full ghosted canvas, Field layout
};

struct RestartData {
    int nx = 0, ny = 0, ng = 0;
    std::vector<RestartField> fields;
    [[nodiscard]] const RestartField* find(std::string_view name) const;
};

using FieldSpan =
    std::span<const std::pair<std::string, const Field<double>*>>;

[[nodiscard]] std::expected<void, Error> write_restart(
    const std::filesystem::path& path, int nx, int ny, int ng, FieldSpan fields);

[[nodiscard]] std::expected<RestartData, Error> read_restart(
    const std::filesystem::path& path);

}  // namespace ns::io
