#include "wavecore/helmholtz/grid.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace {

using wavecore::helmholtz::Grid2D;

int failures = 0;

void check(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void check_close(double actual, double expected, std::string_view message) {
    constexpr double tolerance = 1.0e-14;
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

void test_rectangular_grid() {
    const Grid2D grid{3, 2, 2.0, 3.0};

    check(grid.nx() == 3, "nx must count interior x nodes");
    check(grid.ny() == 2, "ny must count interior y nodes");
    check(grid.unknown_count() == 6, "unknown count must equal nx * ny");
    check_close(grid.lx(), 2.0, "lx must be preserved");
    check_close(grid.ly(), 3.0, "ly must be preserved");
    check_close(grid.hx(), 0.5, "hx must equal lx / (nx + 1)");
    check_close(grid.hy(), 1.0, "hy must equal ly / (ny + 1)");

    check_close(grid.x(0), 0.5, "first interior x coordinate is hx");
    check_close(grid.x(2), 1.5, "last interior x coordinate is nx * hx");
    check_close(grid.y(0), 1.0, "first interior y coordinate is hy");
    check_close(grid.y(1), 2.0, "last interior y coordinate is ny * hy");
}

void test_index_mapping() {
    const Grid2D grid{3, 2, 2.0, 3.0};

    check(grid.linear_index(0, 0) == 0, "first node must map to zero");
    check(grid.linear_index(2, 0) == 2, "x must be the contiguous index");
    check(grid.linear_index(0, 1) == 3, "next y row must start at nx");
    check(grid.linear_index(2, 1) == 5, "last node must map to N - 1");

    for (Grid2D::size_type j = 0; j < grid.ny(); ++j) {
        for (Grid2D::size_type i = 0; i < grid.nx(); ++i) {
            const auto p = grid.linear_index(i, j);
            const auto [mapped_i, mapped_j] = grid.grid_indices(p);
            check(mapped_i == i && mapped_j == j, "index mapping must be bijective");
        }
    }
}

void test_invalid_indices() {
    const Grid2D grid{3, 2, 2.0, 3.0};

    check_throws<std::out_of_range>([&grid] { (void)grid.x(grid.nx()); },
                                    "x(nx) must be rejected");
    check_throws<std::out_of_range>([&grid] { (void)grid.y(grid.ny()); },
                                    "y(ny) must be rejected");
    check_throws<std::out_of_range>([&grid] { (void)grid.linear_index(grid.nx(), 0); },
                                    "i == nx must be rejected");
    check_throws<std::out_of_range>([&grid] { (void)grid.linear_index(0, grid.ny()); },
                                    "j == ny must be rejected");
    check_throws<std::out_of_range>([&grid] { (void)grid.grid_indices(grid.unknown_count()); },
                                    "p == unknown_count must be rejected");
}

void test_invalid_parameters() {
    constexpr double infinity = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    check_throws<std::invalid_argument>([] { (void)Grid2D{0, 1, 1.0, 1.0}; },
                                        "nx == 0 must be rejected");
    check_throws<std::invalid_argument>([] { (void)Grid2D{1, 0, 1.0, 1.0}; },
                                        "ny == 0 must be rejected");
    check_throws<std::invalid_argument>([] { (void)Grid2D{1, 1, 0.0, 1.0}; },
                                        "zero lx must be rejected");
    check_throws<std::invalid_argument>([] { (void)Grid2D{1, 1, 1.0, -1.0}; },
                                        "negative ly must be rejected");
    check_throws<std::invalid_argument>([infinity] { (void)Grid2D{1, 1, infinity, 1.0}; },
                                        "infinite lx must be rejected");
    check_throws<std::invalid_argument>([nan] { (void)Grid2D{1, 1, 1.0, nan}; },
                                        "NaN ly must be rejected");
}

void test_size_overflow() {
    constexpr auto max_size = std::numeric_limits<Grid2D::size_type>::max();
    constexpr auto overflowing_nx = max_size / 2 + 1;

    check_throws<std::overflow_error>(
        [] { (void)Grid2D{overflowing_nx, 2, 1.0, 1.0}; },
        "nx * ny overflow must be rejected");
    check_throws<std::overflow_error>(
        [] { (void)Grid2D{max_size, 1, 1.0, 1.0}; },
        "nx + 1 interval overflow must be rejected");
}

}  // namespace

int main() {
    test_rectangular_grid();
    test_index_mapping();
    test_invalid_indices();
    test_invalid_parameters();
    test_size_overflow();
    return failures == 0 ? 0 : 1;
}
