// cfg2toml: converts a legacy solver configuration pair
// (Input*.cfg + GridTop*.dat) into the canonical TOML schema.
//
//   cfg2toml <legacy.cfg> [out.toml]
//
// With no output argument the TOML is written to stdout. Exit code 0 on
// success, 1 on any parse/validation failure (message on stderr).

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "config/legacy_import.hpp"

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: cfg2toml <legacy.cfg> [out.toml]\n";
        return 1;
    }

    const std::filesystem::path in = argv[1];
    auto cfg = ns::config::import_legacy(in);
    if (!cfg) {
        std::cerr << cfg.error().message;
        return 1;
    }

    const std::string toml_text = ns::config::to_toml(*cfg);

    if (argc == 2) {
        std::cout << toml_text;
        return 0;
    }

    std::ofstream out(argv[2]);
    if (!out) {
        std::cerr << "cannot open output file: " << argv[2] << "\n";
        return 1;
    }
    out << toml_text;
    out.close();
    std::cerr << "wrote " << argv[2] << "\n";
    return 0;
}
