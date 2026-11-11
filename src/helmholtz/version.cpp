#include "wavecore/helmholtz/version.hpp"

#ifndef WAVECORE_VERSION
#error "WAVECORE_VERSION must be defined by the build system"
#endif

namespace wavecore::helmholtz {

std::string_view project_name() noexcept {
    return "WaveCore Research";
}

std::string_view version() noexcept {
    return WAVECORE_VERSION;
}

}  // namespace wavecore::helmholtz
