#include "wavecore/helmholtz/gmres.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {

using namespace wavecore::helmholtz;
using Complex = CsrMatrix::value_type;
using Vector = CsrMatrix::value_array;
using Dense = std::vector<Vector>;
int failures = 0;

void check(bool condition, std::string_view message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

template <typename Callable>
void rejects(Callable callable) {
    try {
        callable();
        check(false, "invalid input accepted");
    } catch (const std::invalid_argument&) {
        // Expected.
    } catch (...) {
        check(false, "unexpected exception type");
    }
}

CsrMatrix sparse(const Dense& dense) {
    CsrMatrix::index_array row_ptr{0};
    CsrMatrix::index_array columns;
    Vector values;
    for (const auto& row : dense) {
        for (std::size_t col = 0; col < row.size(); ++col) {
            if (row[col] != Complex{}) {
                columns.push_back(col);
                values.push_back(row[col]);
            }
        }
        row_ptr.push_back(values.size());
    }
    return {dense.size(), dense.empty() ? 0 : dense.front().size(), row_ptr, columns, values};
}

Vector dense_product(const Dense& a, const Vector& x) {
    Vector b(a.size(), Complex{});
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            b[i] += a[i][j] * x[j];
        }
    }
    return b;
}

double independent_residual(const Dense& a, const Vector& x, const Vector& b) {
    const auto product = dense_product(a, x);
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        numerator += std::norm(b[i] - product[i]);
        denominator += std::norm(b[i]);
    }
    return std::sqrt(numerator) / (denominator == 0.0 ? 1.0 : std::sqrt(denominator));
}

void verify(const Dense& a, const Vector& b, const GmresOptions& options,
            const GmresResult& result, std::string_view label) {
    const double actual = independent_residual(a, result.solution, b);
    check(std::isfinite(actual), "independent residual must be finite");
    check(std::abs(actual - result.final_residual) <=
              64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, actual),
          "reported residual agrees with independent dense computation");
    check(result.iterations <= options.max_iterations, "global iteration limit");
    check(result.residual_history.size() == result.iterations + 1, "true history length");
    check(result.residual_history.back() == result.final_residual, "last history entry");
    for (const auto value : result.residual_history) {
        check(std::isfinite(value) && value >= 0.0, "finite residual history");
    }
    for (const auto value : result.solution) {
        check(std::isfinite(value.real()) && std::isfinite(value.imag()), "finite solution");
    }
    if (result.status == GmresStatus::CONVERGED) {
        check(actual <= options.tolerance, "true residual meets tolerance on success");
    }
    std::cout << label << ": iterations=" << result.iterations
              << ", true_residual=" << std::setprecision(16) << actual
              << ", status=" << static_cast<int>(result.status) << '\n';
}

void solve_known(const Dense& a, const Vector& exact, GmresOptions options,
                 std::string_view label, bool expect_restart = false) {
    const auto b = dense_product(a, exact);
    const Vector initial(exact.size(), Complex{0.125, -0.25});
    const auto result = gmres(sparse(a), b, options, initial);
    verify(a, b, options, result, label);
    check(result.status == GmresStatus::CONVERGED, "known system converges");
    for (std::size_t i = 0; i < exact.size(); ++i) {
        check(std::abs(result.solution[i] - exact[i]) <= 1.0e-8,
              "known solution recovered");
    }
    if (expect_restart) {
        check(result.iterations > options.restart, "test actually crosses restart boundary");
    }
}

void test_systems() {
    solve_known({{{1, 0}, {0, 0}}, {{0, 0}, {1, 0}}},
                {{2, -1}, {-3, 4}}, {1.0e-12, 2, 10}, "identity happy breakdown");
    // Non-Hermitian, genuinely complex, with a dense independent oracle.
    solve_known({{{4, 1}, {1, -2}, {0.5, 0}},
                 {{0, 1}, {5, -1}, {-1, 0.5}},
                 {{2, 0}, {0, -1}, {6, 2}}},
                {{1, 2}, {-2, 0.5}, {0.25, -1}},
                {1.0e-12, 3, 10}, "complex non-Hermitian");
    solve_known({{{2, 1}, {0.2, -0.1}, {0, 0}},
                 {{0, 0}, {3, -0.5}, {0.1, 0.2}},
                 {{0, 0}, {0, 0}, {4, 1}}},
                {{1, -1}, {2, 0.5}, {-1, 2}},
                {1.0e-12, 1, 100}, "complex restart=1", true);
    solve_known({{{2, 0}, {0, 0}, {0, 0}, {0, 0}},
                 {{0, 0}, {3, 0}, {0, 0}, {0, 0}},
                 {{0, 0}, {0, 0}, {4, 0}, {0, 0}},
                 {{0, 0}, {0, 0}, {0, 0}, {5, 0}}},
                {{1, 1}, {2, -1}, {-1, 0.5}, {0.5, 2}},
                {1.0e-12, 2, 100}, "diagonal restart=2", true);
    solve_known({{{0, 0}, {0, 1}}, {{1, 0}, {0, 0}}},
                {{0, 0}, {0, -1}}, {1.0e-12, 2, 5}, "complex permutation");
}

void test_special_cases() {
    const Dense a{{{1, 0}, {0, 0}}, {{0, 0}, {1, 0}}};
    const auto matrix = sparse(a);
    const GmresOptions options{1.0e-12, 2, 10};
    const Vector zero(2, Complex{});
    const Vector identity_rhs{{2, -1}, {-3, 4}};
    auto happy = gmres(matrix, identity_rhs, options);
    verify(a, identity_rhs, options, happy, "identity one-step convergence");
    check(happy.status == GmresStatus::CONVERGED && happy.iterations == 1,
          "identity converges in one step at happy breakdown");

    // With zero x0 and b=e1, the first Arnoldi diagonal is exactly zero.
    const Dense permutation{{{0, 0}, {0, 1}}, {{1, 0}, {0, 0}}};
    const Vector unit_rhs{{1, 0}, {0, 0}};
    auto zero_pivot = gmres(sparse(permutation), unit_rhs, options);
    verify(permutation, unit_rhs, options, zero_pivot, "zero first QR diagonal");
    check(zero_pivot.status == GmresStatus::CONVERGED && zero_pivot.iterations == 2,
          "complex rotation handles zero first diagonal");

    auto result = gmres(matrix, zero, options);
    verify(a, zero, options, result, "zero b and x0");
    check(result.status == GmresStatus::CONVERGED && result.iterations == 0,
          "zero initial residual converges without iteration");
    result = gmres(matrix, zero, options, {{2, 1}, {-1, 3}});
    verify(a, zero, options, result, "zero b nonzero x0");
    check(result.status == GmresStatus::CONVERGED, "absolute zero-rhs stopping rule");
    const Vector exact{{1, 2}, {3, -1}};
    result = gmres(matrix, exact, options, exact);
    verify(a, exact, options, result, "exact x0");
    check(result.status == GmresStatus::CONVERGED && result.iterations == 0,
          "exact x0 needs no iteration");
    const Dense diagonal{{{1, 0}, {0, 0}}, {{0, 0}, {3, 0}}};
    const Vector b{{1, 0}, {1, 0}};
    const GmresOptions limited{1.0e-14, 1, 3};
    result = gmres(sparse(diagonal), b, limited);
    verify(diagonal, b, limited, result, "iteration limit after restarts");
    check(result.status == GmresStatus::NOT_CONVERGED && result.iterations == 3,
          "insufficient global limit returns NOT_CONVERGED");
    const Dense four_diagonal{
        {{1, 0}, {0, 0}, {0, 0}, {0, 0}},
        {{0, 0}, {2, 0}, {0, 0}, {0, 0}},
        {{0, 0}, {0, 0}, {3, 0}, {0, 0}},
        {{0, 0}, {0, 0}, {0, 0}, {4, 0}}};
    const Vector four_rhs(4, Complex{1, 0});
    const GmresOptions partial{1.0e-14, 2, 3};
    result = gmres(sparse(four_diagonal), four_rhs, partial);
    verify(four_diagonal, four_rhs, partial, result, "partial final restart cycle");
    check(result.status == GmresStatus::NOT_CONVERGED && result.iterations == 3,
          "partial restart cycle respects exact global limit");
    const Dense singular{{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}};
    result = gmres(sparse(singular), b, options);
    verify(singular, b, options, result, "inconsistent zero operator");
    check(result.status == GmresStatus::BREAKDOWN, "singular breakdown is not false success");
    result = gmres(CsrMatrix{0, 0, {0}, {}, {}}, {}, options);
    verify({}, {}, options, result, "empty system");
    check(result.status == GmresStatus::CONVERGED && result.iterations == 0,
          "empty system converges");
}

void test_invalid_inputs() {
    const auto matrix = sparse({{{1, 0}, {0, 0}}, {{0, 0}, {1, 0}}});
    const Vector b{{1, 0}, {2, 1}};
    rejects([&] { (void)gmres(CsrMatrix{1, 2, {0, 0}, {}, {}}, {{1, 0}}); });
    rejects([&] { (void)gmres(matrix, {{1, 0}}); });
    rejects([&] { (void)gmres(matrix, b, {}, {{1, 0}}); });
    for (const auto tolerance : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::quiet_NaN()}) {
        rejects([&] { (void)gmres(matrix, b, {tolerance, 2, 10}); });
    }
    rejects([&] { (void)gmres(matrix, b, {1.0e-10, 0, 10}); });
    rejects([&] { (void)gmres(matrix, b, {1.0e-10, 2, 0}); });
    rejects([&] {
        (void)gmres(matrix, {{std::numeric_limits<double>::quiet_NaN(), 0}, {1, 0}});
    });
    rejects([&] {
        (void)gmres(matrix, b, {}, {{1, 0}, {0, std::numeric_limits<double>::infinity()}});
    });
}

}  // namespace

int main() {
    test_systems();
    test_special_cases();
    test_invalid_inputs();
    return failures == 0 ? 0 : 1;
}
