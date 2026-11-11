#include "wavecore/helmholtz/assembly.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using wavecore::helmholtz::CsrMatrix;
using wavecore::helmholtz::Grid2D;
using wavecore::helmholtz::algebraic_manufactured_rhs;
using wavecore::helmholtz::assemble_helmholtz_operator;
using wavecore::helmholtz::continuous_manufactured_rhs;
using wavecore::helmholtz::continuous_mode_eigenvalue;
using wavecore::helmholtz::discrete_mode_eigenvalue;
using wavecore::helmholtz::manufactured_solution;

using Value = CsrMatrix::value_type;
using Vector = CsrMatrix::value_array;

int failures = 0;

void check(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void check_close(Value actual, Value expected, double tolerance, std::string_view message) {
    check(std::abs(actual - expected) <= tolerance, message);
}

void check_close(double actual, double expected, double tolerance, std::string_view message) {
    check(std::abs(actual - expected) <= tolerance, message);
}

template <typename Exception, typename Callable>
void check_throws(Callable&& callable, std::string_view message) {
    try {
        callable();
        check(false, message);
    } catch (const Exception&) {
        // Expected exception.
    } catch (...) {
        check(false, message);
    }
}

std::size_t formula_nnz(std::size_t nx, std::size_t ny) {
    return nx * ny + 2 * (nx - 1) * ny + 2 * nx * (ny - 1);
}

void test_nnz_counts() {
    struct Case {
        std::size_t nx;
        std::size_t ny;
        std::size_t expected;
    };
    constexpr Case cases[]{{1, 1, 1}, {2, 2, 12}, {3, 2, 20}, {5, 5, 105}};

    for (const auto& test_case : cases) {
        const Grid2D grid{test_case.nx, test_case.ny, 1.0, 1.0};
        const auto matrix = assemble_helmholtz_operator(grid, 2.0);
        check(grid.unknown_count() == test_case.nx * test_case.ny,
              "unknown count equals nx * ny");
        check(matrix.nnz() == test_case.expected, "control nnz count");
        check(matrix.nnz() == formula_nnz(test_case.nx, test_case.ny),
              "nnz matches the five-point formula");
    }
}

void test_rectangular_stencil() {
    const Grid2D grid{3, 2, 2.0, 3.0};
    const auto matrix = assemble_helmholtz_operator(grid, 2.0);
    const CsrMatrix::index_array expected_row_ptr{0, 3, 7, 10, 13, 17, 20};
    const CsrMatrix::index_array expected_col_idx{
        0, 1, 3,
        0, 1, 2, 4,
        1, 2, 5,
        0, 3, 4,
        1, 3, 4, 5,
        2, 4, 5};
    const Vector expected_values{
        {6.0, 0.0}, {-4.0, 0.0}, {-1.0, 0.0},
        {-4.0, 0.0}, {6.0, 0.0}, {-4.0, 0.0}, {-1.0, 0.0},
        {-4.0, 0.0}, {6.0, 0.0}, {-1.0, 0.0},
        {-1.0, 0.0}, {6.0, 0.0}, {-4.0, 0.0},
        {-1.0, 0.0}, {-4.0, 0.0}, {6.0, 0.0}, {-4.0, 0.0},
        {-1.0, 0.0}, {-4.0, 0.0}, {6.0, 0.0}};

    check_close(grid.hx(), 0.5, 0.0, "rectangular hx");
    check_close(grid.hy(), 1.0, 0.0, "rectangular hy");
    check(matrix.row_ptr() == expected_row_ptr, "all rectangular row pointers");
    check(matrix.col_idx() == expected_col_idx, "all rectangular column indices");
    check(matrix.values() == expected_values, "all rectangular stencil values");

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        const auto begin = matrix.row_ptr()[row];
        const auto end = matrix.row_ptr()[row + 1];
        check(std::is_sorted(matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(begin),
                             matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(end)),
              "CSR columns are sorted in every row");
        check(std::adjacent_find(
                  matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(begin),
                  matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(end)) ==
                  matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(end),
              "CSR columns are unique in every row");
    }

    const auto row_two_begin = matrix.row_ptr()[2];
    const auto row_two_end = matrix.row_ptr()[3];
    check(std::find(matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(row_two_begin),
                    matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(row_two_end),
                    3) == matrix.col_idx().begin() + static_cast<std::ptrdiff_t>(row_two_end),
          "end of x row does not connect to next x row");
}

void test_symmetry() {
    const Grid2D grid{4, 3, 2.0, 1.5};
    const auto matrix = assemble_helmholtz_operator(grid, 3.0);
    std::vector<Value> dense(matrix.rows() * matrix.cols(), Value{0.0, 0.0});

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t entry = matrix.row_ptr()[row];
             entry < matrix.row_ptr()[row + 1];
             ++entry) {
            dense[row * matrix.cols() + matrix.col_idx()[entry]] = matrix.values()[entry];
        }
    }

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            check(dense[row * matrix.cols() + col] == dense[col * matrix.cols() + row],
                  "assembled matrix equals its transpose");
        }
    }
}

void test_eigenvalues_and_right_hand_sides() {
    constexpr double pi = std::numbers::pi_v<double>;
    const Grid2D one_node{1, 1, 1.0, 1.0};
    check_close(continuous_mode_eigenvalue(one_node, 1, 1),
                2.0 * pi * pi,
                1.0e-13,
                "continuous eigenvalue on the unit square");
    check_close(discrete_mode_eigenvalue(one_node, 1, 1),
                16.0,
                1.0e-13,
                "one-node discrete eigenvalue");

    const Grid2D grid{3, 2, 2.0, 3.0};
    constexpr std::size_t mode_x = 1;
    constexpr std::size_t mode_y = 2;
    constexpr double wavenumber = 0.75;
    const auto exact = manufactured_solution(grid, mode_x, mode_y, Value{1.0, -0.5});
    const auto continuous_rhs = continuous_manufactured_rhs(
        grid, mode_x, mode_y, wavenumber, Value{1.0, -0.5});
    const auto algebraic_rhs = algebraic_manufactured_rhs(
        grid, mode_x, mode_y, wavenumber, Value{1.0, -0.5});
    const double continuous_scale =
        continuous_mode_eigenvalue(grid, mode_x, mode_y) - wavenumber * wavenumber;
    const double discrete_scale =
        discrete_mode_eigenvalue(grid, mode_x, mode_y) - wavenumber * wavenumber;

    check(std::abs(continuous_scale - discrete_scale) > 1.0e-6,
          "continuous and discrete eigenvalues differ on a finite grid");
    for (std::size_t p = 0; p < exact.size(); ++p) {
        check_close(continuous_rhs[p], continuous_scale * exact[p], 1.0e-13,
                    "continuous manufactured right-hand side");
        check_close(algebraic_rhs[p], discrete_scale * exact[p], 1.0e-13,
                    "algebraic manufactured right-hand side");
    }
}

void check_algebraic_identity(const Grid2D& grid,
                              std::size_t mode_x,
                              std::size_t mode_y,
                              double wavenumber) {
    const auto matrix = assemble_helmholtz_operator(grid, wavenumber);
    const auto exact = manufactured_solution(grid, mode_x, mode_y, Value{1.0, 0.25});
    const auto expected =
        algebraic_manufactured_rhs(grid, mode_x, mode_y, wavenumber, Value{1.0, 0.25});
    const auto actual = matrix.multiply(exact);

    const double operator_scale =
        4.0 / (grid.hx() * grid.hx()) +
        4.0 / (grid.hy() * grid.hy()) + wavenumber * wavenumber;
    const double tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * operator_scale;

    for (std::size_t p = 0; p < actual.size(); ++p) {
        check_close(actual[p], expected[p], tolerance,
                    "A_h times exact mode equals algebraic right-hand side");
    }
}

void test_algebraic_identity() {
    check_algebraic_identity(Grid2D{3, 2, 2.0, 3.0}, 1, 2, 0.75);
    check_algebraic_identity(Grid2D{5, 4, 1.0, 2.0}, 2, 3, 4.0);
    check_algebraic_identity(Grid2D{7, 6, 3.0, 1.5}, 3, 2, 8.0);
}

void test_invalid_parameters_and_overflow() {
    const Grid2D grid{3, 2, 2.0, 3.0};
    constexpr double infinity = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    check_throws<std::invalid_argument>(
        [&grid] { (void)assemble_helmholtz_operator(grid, -1.0); },
        "negative wavenumber must be rejected");
    check_throws<std::invalid_argument>(
        [&grid, infinity] { (void)assemble_helmholtz_operator(grid, infinity); },
        "infinite wavenumber must be rejected");
    check_throws<std::invalid_argument>(
        [&grid, nan] { (void)assemble_helmholtz_operator(grid, nan); },
        "NaN wavenumber must be rejected");
    check_throws<std::overflow_error>(
        [&grid] {
            (void)assemble_helmholtz_operator(
                grid, std::sqrt(std::numeric_limits<double>::max()) * 2.0);
        },
        "squared wavenumber overflow must be rejected");
    check_throws<std::invalid_argument>(
        [&grid] { (void)manufactured_solution(grid, 0, 1); },
        "zero mode number must be rejected");
    check_throws<std::invalid_argument>(
        [&grid] { (void)manufactured_solution(grid, 4, 1); },
        "mode beyond nx must be rejected");
    check_throws<std::invalid_argument>(
        [&grid, infinity] { (void)manufactured_solution(grid, 1, 1, {infinity, 0.0}); },
        "non-finite amplitude must be rejected");

    const Grid2D tiny_spacing{1, 1, std::numeric_limits<double>::min(), 1.0};
    check_throws<std::overflow_error>(
        [&tiny_spacing] { (void)assemble_helmholtz_operator(tiny_spacing, 1.0); },
        "non-finite inverse squared spacing must be rejected");

    constexpr auto max_size = std::numeric_limits<std::size_t>::max();
    const Grid2D huge_grid{max_size / 4, 2, 1.0, 1.0};
    check_throws<std::overflow_error>(
        [&huge_grid] { (void)assemble_helmholtz_operator(huge_grid, 1.0); },
        "CSR nnz overflow must be rejected before allocation");
}

}  // namespace

int main() {
    test_nnz_counts();
    test_rectangular_stencil();
    test_symmetry();
    test_eigenvalues_and_right_hand_sides();
    test_algebraic_identity();
    test_invalid_parameters_and_overflow();
    return failures == 0 ? 0 : 1;
}
