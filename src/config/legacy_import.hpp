#pragma once

// Imports a legacy solver configuration (Input*.cfg + its GridTop*.dat
// topology companion) into the typed ns::config::Config model.
//
// Topology resolution: the GridTop file is looked up next to the .cfg file
// (the legacy code hard-coded "../input/"; the importer uses the config's own
// directory instead). The legacy 'bound_cell' column is ignored -- it encodes
// ghost-row indices that are derivable from the edge in the new model.

#include <expected>
#include <filesystem>

#include "config/config.hpp"

namespace ns::config {

std::expected<Config, Error> import_legacy(const std::filesystem::path& legacy_cfg);

}  // namespace ns::config
