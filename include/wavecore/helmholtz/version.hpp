#pragma once

#include <string_view>

namespace wavecore::helmholtz {

[[nodiscard]] std::string_view project_name() noexcept;
[[nodiscard]] std::string_view version() noexcept;

}  // namespace wavecore::helmholtz
