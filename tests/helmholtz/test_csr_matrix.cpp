#include "wavecore/helmholtz/csr_matrix.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using Matrix = wavecore::helmholtz::CsrMatrix;
using Value = Matrix::value_type;
using Vector = Matrix::value_array;

int failures = 0;

void check(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

void check_close(Value actual, Value expected, std::string_view message) {
    constexpr double tolerance = 1.0e-13;
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

Matrix make_known_matrix() {
    return Matrix{
        2,
        2,
        {0, 2, 4},
        {0, 1, 0, 1},
        {Value{1.0, 2.0}, Value{3.0, -1.0}, Value{-2.0, 0.5}, Value{4.0, 0.0}}};
}

void test_known_complex_product() {
    const auto matrix = make_known_matrix();
    const Vector input{Value{2.0, -1.0}, Value{-1.0, 3.0}};
    const auto result = matrix.multiply(input);

    check(matrix.rows() == 2, "known matrix row count");
    check(matrix.cols() == 2, "known matrix column count");
    check(matrix.nnz() == 4, "known matrix nnz");
    check(result.size() == 2, "known product size");
    check_close(result[0], Value{4.0, 13.0}, "known product first entry");
    check_close(result[1], Value{-7.5, 15.0}, "known product second entry");
}

void test_rectangular_matrix() {
    const Matrix matrix{
        2,
        3,
        {0, 2, 3},
        {1, 2, 0},
        {Value{0.0, 2.0}, Value{1.0, 0.0}, Value{-1.0, 0.0}}};
    const Vector input{Value{1.0, 1.0}, Value{2.0, -1.0}, Value{-3.0, 2.0}};
    const auto result = matrix.multiply(input);

    check(result.size() == 2, "rectangular product size equals rows");
    check_close(result[0], Value{-1.0, 6.0}, "rectangular product first entry");
    check_close(result[1], Value{-1.0, -1.0}, "rectangular product second entry");
}

void test_empty_matrix_and_rows() {
    const Matrix empty{0, 0, {0}, {}, {}};
    check(empty.rows() == 0 && empty.cols() == 0 && empty.nnz() == 0,
          "zero-by-zero matrix metadata");
    check(empty.multiply({}).empty(), "zero-by-zero product must be empty");

    const Matrix empty_rows{3, 2, {0, 0, 0, 0}, {}, {}};
    const Vector input{Value{1.0, 2.0}, Value{-3.0, 4.0}};
    const auto result = empty_rows.multiply(input);
    check(result == Vector(3, Value{0.0, 0.0}), "empty rows must produce zero entries");

    const Matrix one_empty_row{2, 2, {0, 1, 1}, {1}, {Value{2.0, 0.0}}};
    const auto mixed_result = one_empty_row.multiply(input);
    check_close(mixed_result[0], Value{-6.0, 8.0}, "nonempty row around an empty row");
    check_close(mixed_result[1], Value{0.0, 0.0}, "explicitly empty row result");
}

void test_invalid_structures() {
    constexpr double infinity = std::numeric_limits<double>::infinity();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 1, {0}, {}, {}}; }, "row_ptr size mismatch");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 1, {1, 1}, {0}, {Value{1.0, 0.0}}}; },
        "row_ptr must start at zero");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 1, {0, 0}, {0}, {Value{1.0, 0.0}}}; },
        "row_ptr must end at nnz");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 1, {0, 1}, {}, {Value{1.0, 0.0}}}; },
        "column and value array sizes must match");
    check_throws<std::invalid_argument>(
        [] {
            (void)Matrix{3,
                         2,
                         {0, 1, 0, 2},
                         {0, 1},
                         {Value{1.0, 0.0}, Value{2.0, 0.0}}};
        },
        "row_ptr must not decrease");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 2, {0, 2}, {1, 1}, {Value{1.0, 0.0}, Value{2.0, 0.0}}}; },
        "duplicate column indices must be rejected");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 2, {0, 2}, {1, 0}, {Value{1.0, 0.0}, Value{2.0, 0.0}}}; },
        "descending column indices must be rejected");
    check_throws<std::invalid_argument>(
        [] { (void)Matrix{1, 2, {0, 1}, {2}, {Value{1.0, 0.0}}}; },
        "out-of-range column index must be rejected");
    check_throws<std::invalid_argument>(
        [infinity] { (void)Matrix{1, 1, {0, 1}, {0}, {Value{infinity, 0.0}}}; },
        "infinite real value must be rejected");
    check_throws<std::invalid_argument>(
        [nan] { (void)Matrix{1, 1, {0, 1}, {0}, {Value{0.0, nan}}}; },
        "NaN imaginary value must be rejected");
}

void test_size_and_vector_errors() {
    constexpr auto max_size = std::numeric_limits<Matrix::size_type>::max();

    check_throws<std::overflow_error>(
        [] { (void)Matrix{max_size, 0, {0}, {}, {}}; },
        "rows + 1 overflow must be rejected");

    const auto matrix = make_known_matrix();
    check_throws<std::invalid_argument>(
        [&matrix] { (void)matrix.multiply(Vector{Value{1.0, 0.0}}); },
        "incorrect input vector size must be rejected");
}

void test_multiplication_preserves_inputs() {
    const auto matrix = make_known_matrix();
    const auto original_row_ptr = matrix.row_ptr();
    const auto original_col_idx = matrix.col_idx();
    const auto original_values = matrix.values();
    Vector input{Value{2.0, -1.0}, Value{-1.0, 3.0}};
    const auto original_input = input;

    (void)matrix.multiply(input);

    check(matrix.row_ptr() == original_row_ptr, "multiply must preserve row_ptr");
    check(matrix.col_idx() == original_col_idx, "multiply must preserve col_idx");
    check(matrix.values() == original_values, "multiply must preserve values");
    check(input == original_input, "multiply must preserve the input vector");
}

}  // namespace

int main() {
    test_known_complex_product();
    test_rectangular_matrix();
    test_empty_matrix_and_rows();
    test_invalid_structures();
    test_size_and_vector_errors();
    test_multiplication_preserves_inputs();
    return failures == 0 ? 0 : 1;
}
