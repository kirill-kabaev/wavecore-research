#include "wavecore/helmholtz/gmres.hpp"
#include "../../src/helmholtz/gmres_diagnostics.hpp"

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

// Independent dense arithmetic in wider precision avoids squaring tiny/huge
// doubles and reduces cancellation in the verification oracle.
double wide_residual(const Dense& a, const Vector& x, const Vector& b) {
    using Wide = std::complex<long double>;
    long double error = 0.0L;
    long double rhs_norm = 0.0L;
    for (std::size_t i = 0; i < b.size(); ++i) {
        Wide product{};
        for (std::size_t j = 0; j < x.size(); ++j) {
            product += Wide(a[i][j]) * Wide(x[j]);
        }
        error = std::hypot(error, std::abs(Wide(b[i]) - product));
        rhs_norm = std::hypot(rhs_norm, std::abs(Wide(b[i])));
    }
    return static_cast<double>(error / (rhs_norm == 0.0L ? 1.0L : rhs_norm));
}

double relative_error(const Vector& x, const Vector& exact) {
    long double error = 0.0L;
    long double reference = 0.0L;
    for (std::size_t i = 0; i < x.size(); ++i) {
        error = std::hypot(error, std::abs(std::complex<long double>(x[i]) -
                                         std::complex<long double>(exact[i])));
        reference = std::hypot(reference, std::abs(std::complex<long double>(exact[i])));
    }
    return static_cast<double>(error / reference);
}

struct Diagnostic {
    double max_defect = 0.0;
    double last_estimate = 0.0;
    double max_residual_gap = 0.0;
    std::size_t samples = 0;
    std::size_t max_basis_size = 0;
    bool optimistic_estimate = false;

    detail::ArnoldiObserver observer(double tolerance) {
        return [this, tolerance](const std::vector<Vector>& basis, const Vector& next,
                                  double next_norm, bool breakdown,
                                  double estimated, double actual) {
            auto vectors = basis;
            if (!breakdown) {
                auto normalized = next;
                for (auto& entry : normalized) {
                    entry /= next_norm;
                }
                vectors.push_back(normalized);
            }
            // Frobenius norm ||V^H V - I||, independently accumulated in long double.
            long double defect = 0.0L;
            for (std::size_t i = 0; i < vectors.size(); ++i) {
                for (std::size_t j = 0; j < vectors.size(); ++j) {
                    std::complex<long double> product{};
                    for (std::size_t k = 0; k < next.size(); ++k) {
                        product += std::conj(std::complex<long double>(vectors[i][k])) *
                                   std::complex<long double>(vectors[j][k]);
                    }
                    if (i == j) {
                        product -= 1.0L;
                    }
                    defect = std::hypot(defect, std::abs(product));
                }
            }
            max_defect = std::max(max_defect, static_cast<double>(defect));
            max_basis_size = std::max(max_basis_size, vectors.size());
            last_estimate = estimated;
            max_residual_gap = std::max(max_residual_gap, std::abs(estimated - actual));
            optimistic_estimate |= estimated <= tolerance && actual > tolerance;
            ++samples;
            check(std::isfinite(estimated) && std::isfinite(actual), "finite diagnostic residuals");
        };
    }
};

GmresResult diagnostic_solve(const Dense& a, const Vector& exact,
                              const GmresOptions& options, std::string_view label,
                              Diagnostic& diagnostic) {
    const auto b = dense_product(a, exact);
    const auto result = detail::gmres_observed(sparse(a), b, options, {},
                                              diagnostic.observer(options.tolerance));
    const auto public_result = gmres(sparse(a), b, options);
    check(public_result.status == result.status &&
              public_result.iterations == result.iterations &&
              public_result.solution == result.solution &&
              public_result.residual_history == result.residual_history,
          "internal observation leaves public solver behavior unchanged");
    const double actual = wide_residual(a, result.solution, b);
    const double error = relative_error(result.solution, exact);
    check(std::isfinite(actual) && std::isfinite(error), "finite independent error and residual");
    check(result.iterations <= options.max_iterations, "diagnostic global iteration limit");
    check(result.residual_history.size() == result.iterations + 1, "diagnostic history size");
    for (const auto entry : result.solution) {
        check(std::isfinite(entry.real()) && std::isfinite(entry.imag()), "finite diagnostic solution");
    }
    if (result.status == GmresStatus::CONVERGED) {
        check(actual <= options.tolerance, "independent wide residual confirms success");
    }
    const double comparison_tolerance = 128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, actual);
    check(std::abs(actual - result.final_residual) <= comparison_tolerance,
          "wide and reported residual agree within roundoff");
    check(diagnostic.samples > 0, "actual Arnoldi basis observed");
    check(diagnostic.max_defect < 1.0e-12, "actual two-pass MGS basis remains orthogonal");
    std::cout << label << ": iterations=" << result.iterations
              << ", status=" << static_cast<int>(result.status)
              << ", error=" << error << ", wide_residual=" << actual
              << ", estimated=" << diagnostic.last_estimate
              << ", orthogonality_F=" << diagnostic.max_defect
              << ", max_residual_gap=" << diagnostic.max_residual_gap << '\n';
    return result;
}

double single_pass_arnoldi_defect(const Dense& a, const Vector& b, std::size_t steps) {
    // Diagnostic control only: a local one-pass implementation deliberately omits
    // the production solver's reorthogonalization, QR and stopping rules.
    auto vector_norm = [](const Vector& values) {
        double result = 0.0;
        for (auto value : values) {
            result = std::hypot(result, std::abs(value));
        }
        return result;
    };
    Vector initial = b;
    const double beta = vector_norm(initial);
    for (auto& value : initial) {
        value /= beta;
    }
    std::vector<Vector> basis{initial};
    for (std::size_t column = 0; column < steps; ++column) {
        auto w = dense_product(a, basis[column]);
        for (const auto& v : basis) {
            Complex projection{};
            for (std::size_t i = 0; i < w.size(); ++i) {
                projection += std::conj(v[i]) * w[i];
            }
            for (std::size_t i = 0; i < w.size(); ++i) {
                w[i] -= projection * v[i];
            }
        }
        const double length = vector_norm(w);
        check(length > 0.0 && std::isfinite(length), "single-pass diagnostic has a direction");
        if (length == 0.0 || !std::isfinite(length)) {
            break;
        }
        for (auto& value : w) {
            value /= length;
        }
        basis.push_back(w);
    }
    long double defect = 0.0L;
    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = 0; j < basis.size(); ++j) {
            std::complex<long double> product{};
            for (std::size_t k = 0; k < b.size(); ++k) {
                product += std::conj(std::complex<long double>(basis[i][k])) *
                           std::complex<long double>(basis[j][k]);
            }
            product -= (i == j ? 1.0L : 0.0L);
            defect = std::hypot(defect, std::abs(product));
        }
    }
    return static_cast<double>(defect);
}

void test_extended_systems() {
    // Unitary diagonal similarity of Hilbert(10): complex Hermitian, ill-conditioned.
    const std::size_t n = 10;
    Dense hilbert(n, Vector(n));
    Vector exact(n);
    for (std::size_t i = 0; i < n; ++i) {
        exact[i] = {1.0, static_cast<double>(i % 3) * 0.25};
        for (std::size_t j = 0; j < n; ++j) {
            const double angle = 0.3 * (static_cast<double>(i) - static_cast<double>(j));
            hilbert[i][j] = std::polar(1.0 / static_cast<double>(i + j + 1), angle);
        }
    }
    Diagnostic hilbert_diagnostic;
    const auto hilbert_result = diagnostic_solve(hilbert, exact, {1.0e-12, n, 30},
                                                 "complex Hilbert(10)", hilbert_diagnostic);
    check(hilbert_result.status == GmresStatus::CONVERGED, "Hilbert residual converges");
    check(hilbert_diagnostic.max_basis_size >= 7, "Hilbert exercises a long Arnoldi basis");
    const double single_pass_defect = single_pass_arnoldi_defect(
        hilbert, dense_product(hilbert, exact), 8);
    std::cout << "Hilbert single-pass MGS control: orthogonality_F=" << single_pass_defect << '\n';
    check(single_pass_defect > 1.0e-8,
          "diagnostic system exposes loss of orthogonality without reorthogonalization");

    // Exact condition number is 1e14. Small residual need not resolve the weak mode.
    const Dense near_singular{{{1, 0}, {0, 0}}, {{0, 0}, {1.0e-14, 0}}};
    Diagnostic near_diagnostic;
    const Vector near_exact{{1, 0}, {1, 0}};
    const auto near_result = diagnostic_solve(near_singular, near_exact, {1.0e-12, 2, 4},
                                              "near-singular diagonal", near_diagnostic);
    check(near_result.status == GmresStatus::CONVERGED, "weak mode below residual tolerance");
    check(relative_error(near_result.solution, near_exact) > 0.5,
          "small residual demonstrably does not guarantee solution accuracy");

    // A fixed complex nonnormal Jordan-like matrix, with independently prescribed x.
    Dense jordan(6, Vector(6, Complex{}));
    Vector jordan_exact(6);
    for (std::size_t i = 0; i < 6; ++i) {
        jordan[i][i] = {1, 0.5};
        if (i + 1 < 6) {
            jordan[i][i + 1] = {2, -1};
        }
        jordan_exact[i] = {1.0 + static_cast<double>(i), -0.25};
    }
    Diagnostic jordan_diagnostic;
    const auto jordan_result = diagnostic_solve(jordan, jordan_exact, {1.0e-12, 6, 12},
                                                "complex nonnormal Jordan", jordan_diagnostic);
    check(jordan_result.status == GmresStatus::CONVERGED, "nonnormal system converges");
    check(relative_error(jordan_result.solution, jordan_exact) < 1.0e-10,
          "nonnormal known solution recovered");

    const Dense base{{{4, 1}, {1, -2}, {0.5, 0}},
                     {{0, 1}, {5, -1}, {-1, 0.5}},
                     {{2, 0}, {0, -1}, {6, 2}}};
    const Vector base_exact{{1, 2}, {-2, 0.5}, {0.25, -1}};
    for (const double scale : {1.0e-200, 1.0, 1.0e200}) {
        auto scaled = base;
        for (auto& row : scaled) {
            for (auto& entry : row) {
                entry *= scale;
            }
        }
        Diagnostic diagnostic;
        std::cout << "scale=" << scale << ' ';
        const auto result = diagnostic_solve(scaled, base_exact, {1.0e-12, 99, 6},
                                             "scaled complex system", diagnostic);
        check(result.status == GmresStatus::CONVERGED, "scaled system converges");
        check(relative_error(result.solution, base_exact) < 1.0e-10,
              "scaled solution accuracy");
        check(result.iterations <= 3, "restart >= n capped at n");
    }

    Diagnostic gap_diagnostic;
    const auto gap_result = diagnostic_solve(base, base_exact, {1.0e-30, 3, 6},
                                             "estimated/true residual gap", gap_diagnostic);
    check(gap_diagnostic.optimistic_estimate, "estimated convergence with unmet true tolerance observed");
    check(gap_result.status != GmresStatus::CONVERGED, "optimistic estimate cannot trigger success");
}

void test_stagnation_and_limits() {
    // Invertible skew-Hermitian A, b=e1: GMRES(1) stagnates exactly at x=0.
    const Dense skew{{{0, 0}, {-1, 0}}, {{1, 0}, {0, 0}}};
    const Vector exact{{0, 0}, {-1, 0}};
    Diagnostic diagnostic;
    const auto result = diagnostic_solve(skew, exact, {1.0e-12, 1, 5},
                                         "exact GMRES(1) stagnation", diagnostic);
    check(result.status == GmresStatus::NOT_CONVERGED && result.iterations == 5,
          "stagnation reaches limit without false success");
    for (const auto residual : result.residual_history) {
        check(residual == 1.0, "stagnation history remains exactly one");
    }
    for (const auto limit : {std::size_t{1}, std::size_t{2}}) {
        Diagnostic boundary;
        const auto boundary_result = diagnostic_solve(skew, exact, {1.0e-12, 2, limit},
                                                       "iteration boundary", boundary);
        check(boundary_result.iterations == limit, "exact boundary iteration count");
        check(boundary_result.status == (limit == 1 ? GmresStatus::NOT_CONVERGED :
                                                     GmresStatus::CONVERGED),
              "success exactly at limit, failure before it");
    }
    bool overflow = false;
    try {
        // Finite A and x0, but their product cannot be represented as double.
        const double largest = std::numeric_limits<double>::max();
        (void)gmres(sparse({{{largest, 0}}}), {{1, 0}}, {}, {{2, 0}});
    } catch (const std::overflow_error&) {
        overflow = true;
    }
    check(overflow, "nonfinite arithmetic raises an error instead of returning success");
}

}  // namespace

int main() {
    test_systems();
    test_special_cases();
    test_invalid_inputs();
    test_extended_systems();
    test_stagnation_and_limits();
    return failures == 0 ? 0 : 1;
}
