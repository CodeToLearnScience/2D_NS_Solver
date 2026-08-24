#pragma once

// Shared lightweight error type for the new solver infrastructure.
// (config::Error is kept separate for now; consolidation planned when the
// config layer is next touched.)

#include <format>
#include <string>

namespace ns {

struct Error {
    std::string message;
};

}  // namespace ns
