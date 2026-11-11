#include "wavecore/helmholtz/assembly.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace wavecore::helmholtz {
namespace {

using size_type = CsrMatrix::size_type;
using value_type = CsrMatrix::value_type;
using value_array = CsrMatrix::value_array;

[[nodiscard]] size_type checked_add(size_type left, size_type right) {
    constexpr auto max_size = std::numeric_limits<size_type>::max();
    if (left > max_size - right) {
        throw std::overflow_error("Sparse matrix size addition overflow");
    }
    return left + right;
}

[[nodiscard]] size_type checked_multiply(size_type left, size_type right) {
    constexpr auto max_size = std::numeric_limits<size_type>::max();
    if (left != 0 && right > max_size / left) {
        throw std::overflow_error("Sparse matrix size multiplication overflow");
    }
    return left * right;
}

[[nodiscard]] size_type expected_nnz(const Grid2D& grid) {
    const auto horizontal_edges = checked_multiply(grid.nx() - 1, grid.ny());
    const auto vertical_edges = checked_multiply(grid.nx(), grid.ny() - 1);
    const auto off_diagonal = checked_multiply(
        2, checked_add(horizontal_edges, vertical_edges));
    return checked_add(grid.unknown_count(), off_diagonal);
}

void validate_wavenumber(double wavenumber) {
    if (!std::isfinite(wavenumber) || wavenumber < 0.0) {
        throw std::invalid_argument("Wavenumber must be finite and nonnegative");
    }
}

[[nodiscard]] double wavenumber_squared(double wavenumber) {
    validate_wavenumber(wavenumber);
    const double result = wavenumber * wavenumber;
    if (!std::isfinite(result)) {
        throw std::overflow_error("Squared wavenumber is not finite");
    }
    return result;
}

void validate_mode(const Grid2D& grid, size_type mode_x, size_type mode_y) {
    if (mode_x == 0 || mode_y == 0) {
        throw std::invalid_argument("Manufactured mode numbers must be positive");
    }
    if (mode_x > grid.nx() || mode_y > grid.ny()) {
        throw std::invalid_argument("Manufactured mode numbers must fit the interior grid");
    }
}

void validate_amplitude(value_type amplitude) {
    if (!std::isfinite(amplitude.real()) || !std::isfinite(amplitude.imag())) {
        throw std::invalid_argument("Manufactured solution amplitude must be finite");
    }
}

[[nodiscard]] double inverse_spacing_squared(double spacing) {
    const double inverse_spacing = 1.0 / spacing;
    const double result = inverse_spacing * inverse_spacing;
    if (!std::isfinite(result)) {
        throw std::overflow_error("Inverse squared grid spacing is not finite");
    }
    return result;
}

[[nodiscard]] value_array scaled_solution(const Grid2D& grid,
                                          size_type mode_x,
                                          size_type mode_y,
                                          double scale,
                                          value_type amplitude) {
    auto result = manufactured_solution(grid, mode_x, mode_y, amplitude);
    for (auto& value : result) {
        value *= scale;
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::overflow_error("Manufactured right-hand side is not finite");
        }
    }
    return result;
}

}  // namespace

CsrMatrix assemble_helmholtz_operator(const Grid2D& grid, double wavenumber) {
    const double inverse_hx_squared = inverse_spacing_squared(grid.hx());
    const double inverse_hy_squared = inverse_spacing_squared(grid.hy());
    const double k_squared = wavenumber_squared(wavenumber);

    const double x_neighbor = -inverse_hx_squared;
    const double y_neighbor = -inverse_hy_squared;
    const double diagonal =
        2.0 * inverse_hx_squared + 2.0 * inverse_hy_squared - k_squared;
    if (!std::isfinite(x_neighbor) || !std::isfinite(y_neighbor) ||
        !std::isfinite(diagonal)) {
        throw std::overflow_error("Helmholtz stencil coefficient is not finite");
    }

    const auto nonzero_count = expected_nnz(grid);
    const auto row_ptr_size = checked_add(grid.unknown_count(), 1);

    CsrMatrix::index_array row_ptr;
    CsrMatrix::index_array col_idx;
    CsrMatrix::value_array values;
    if (row_ptr_size > row_ptr.max_size() || nonzero_count > col_idx.max_size() ||
        nonzero_count > values.max_size()) {
        throw std::overflow_error("Helmholtz CSR structure exceeds container limits");
    }

    row_ptr.reserve(row_ptr_size);
    col_idx.reserve(nonzero_count);
    values.reserve(nonzero_count);
    row_ptr.push_back(0);

    for (size_type j = 0; j < grid.ny(); ++j) {
        for (size_type i = 0; i < grid.nx(); ++i) {
            const auto p = grid.linear_index(i, j);

            if (j > 0) {
                col_idx.push_back(p - grid.nx());
                values.emplace_back(y_neighbor, 0.0);
            }
            if (i > 0) {
                col_idx.push_back(p - 1);
                values.emplace_back(x_neighbor, 0.0);
            }

            col_idx.push_back(p);
            values.emplace_back(diagonal, 0.0);

            if (i + 1 < grid.nx()) {
                col_idx.push_back(p + 1);
                values.emplace_back(x_neighbor, 0.0);
            }
            if (j + 1 < grid.ny()) {
                col_idx.push_back(p + grid.nx());
                values.emplace_back(y_neighbor, 0.0);
            }

            row_ptr.push_back(col_idx.size());
        }
    }

    if (col_idx.size() != nonzero_count || values.size() != nonzero_count) {
        throw std::logic_error("Assembled Helmholtz matrix has an unexpected nnz count");
    }

    return CsrMatrix{grid.unknown_count(),
                     grid.unknown_count(),
                     std::move(row_ptr),
                     std::move(col_idx),
                     std::move(values)};
}

double continuous_mode_eigenvalue(const Grid2D& grid,
                                  size_type mode_x,
                                  size_type mode_y) {
    validate_mode(grid, mode_x, mode_y);
    constexpr double pi = std::numbers::pi_v<double>;
    const double x_wavenumber = static_cast<double>(mode_x) * pi / grid.lx();
    const double y_wavenumber = static_cast<double>(mode_y) * pi / grid.ly();
    const double result =
        x_wavenumber * x_wavenumber + y_wavenumber * y_wavenumber;
    if (!std::isfinite(result)) {
        throw std::overflow_error("Continuous mode eigenvalue is not finite");
    }
    return result;
}

double discrete_mode_eigenvalue(const Grid2D& grid,
                                size_type mode_x,
                                size_type mode_y) {
    validate_mode(grid, mode_x, mode_y);
    constexpr double pi = std::numbers::pi_v<double>;
    const double x_angle =
        static_cast<double>(mode_x) * pi /
        (2.0 * (static_cast<double>(grid.nx()) + 1.0));
    const double y_angle =
        static_cast<double>(mode_y) * pi /
        (2.0 * (static_cast<double>(grid.ny()) + 1.0));
    const double x_sine = std::sin(x_angle);
    const double y_sine = std::sin(y_angle);
    const double result =
        4.0 * inverse_spacing_squared(grid.hx()) * x_sine * x_sine +
        4.0 * inverse_spacing_squared(grid.hy()) * y_sine * y_sine;
    if (!std::isfinite(result)) {
        throw std::overflow_error("Discrete mode eigenvalue is not finite");
    }
    return result;
}

value_array manufactured_solution(const Grid2D& grid,
                                  size_type mode_x,
                                  size_type mode_y,
                                  value_type amplitude) {
    validate_mode(grid, mode_x, mode_y);
    validate_amplitude(amplitude);
    constexpr double pi = std::numbers::pi_v<double>;

    value_array result(grid.unknown_count());
    for (size_type j = 0; j < grid.ny(); ++j) {
        const double y_phase =
            (static_cast<double>(j) + 1.0) * static_cast<double>(mode_y) * pi /
            (static_cast<double>(grid.ny()) + 1.0);
        const double y_factor = std::sin(y_phase);

        for (size_type i = 0; i < grid.nx(); ++i) {
            const double x_phase =
                (static_cast<double>(i) + 1.0) * static_cast<double>(mode_x) * pi /
                (static_cast<double>(grid.nx()) + 1.0);
            result[grid.linear_index(i, j)] = amplitude * (std::sin(x_phase) * y_factor);
        }
    }
    return result;
}

value_array continuous_manufactured_rhs(const Grid2D& grid,
                                        size_type mode_x,
                                        size_type mode_y,
                                        double wavenumber,
                                        value_type amplitude) {
    const double scale = continuous_mode_eigenvalue(grid, mode_x, mode_y) -
                         wavenumber_squared(wavenumber);
    return scaled_solution(grid, mode_x, mode_y, scale, amplitude);
}

value_array algebraic_manufactured_rhs(const Grid2D& grid,
                                       size_type mode_x,
                                       size_type mode_y,
                                       double wavenumber,
                                       value_type amplitude) {
    const double scale = discrete_mode_eigenvalue(grid, mode_x, mode_y) -
                         wavenumber_squared(wavenumber);
    return scaled_solution(grid, mode_x, mode_y, scale, amplitude);
}

}  // namespace wavecore::helmholtz
