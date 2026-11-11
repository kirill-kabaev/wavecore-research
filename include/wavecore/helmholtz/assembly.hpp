#pragma once

#include "wavecore/helmholtz/csr_matrix.hpp"
#include "wavecore/helmholtz/grid.hpp"

#include <cstddef>

namespace wavecore::helmholtz {

[[nodiscard]] CsrMatrix assemble_helmholtz_operator(const Grid2D& grid,
                                                    double wavenumber);

[[nodiscard]] double continuous_mode_eigenvalue(const Grid2D& grid,
                                                std::size_t mode_x,
                                                std::size_t mode_y);

[[nodiscard]] double discrete_mode_eigenvalue(const Grid2D& grid,
                                              std::size_t mode_x,
                                              std::size_t mode_y);

[[nodiscard]] CsrMatrix::value_array manufactured_solution(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    CsrMatrix::value_type amplitude = {1.0, 0.0});

[[nodiscard]] CsrMatrix::value_array continuous_manufactured_rhs(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    double wavenumber,
    CsrMatrix::value_type amplitude = {1.0, 0.0});

[[nodiscard]] CsrMatrix::value_array algebraic_manufactured_rhs(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    double wavenumber,
    CsrMatrix::value_type amplitude = {1.0, 0.0});

}  // namespace wavecore::helmholtz
