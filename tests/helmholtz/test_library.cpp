#include "wavecore/helmholtz/version.hpp"

#include <iostream>
#include <string_view>

#ifndef WAVECORE_EXPECTED_VERSION
#error "WAVECORE_EXPECTED_VERSION must be defined by the build system"
#endif

int main() {
    constexpr std::string_view expected_name{"WaveCore Research"};
    constexpr std::string_view expected_version{WAVECORE_EXPECTED_VERSION};

    if (wavecore::helmholtz::project_name() != expected_name) {
        std::cerr << "Unexpected project name\n";
        return 1;
    }

    if (wavecore::helmholtz::version() != expected_version) {
        std::cerr << "Unexpected project version\n";
        return 1;
    }

    return 0;
}
