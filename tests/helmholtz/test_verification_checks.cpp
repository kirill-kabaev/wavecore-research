#include "verification_checks.hpp"
#include "wavecore/helmholtz/assembly.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

namespace checks = wavecore::test::verification;
using checks::Profile;
using checks::Reason;
using namespace wavecore::helmholtz;

std::size_t passed = 0;
constexpr GmresOptions options{1.0e-11, 20, 200};
constexpr long double fraction = 0.001L;
const long double order_slack = std::log2((1.0L + fraction) / (1.0L - fraction));

void test_assert(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

template<class Callable>
void accepted(std::string_view label, Callable&& callable) {
    callable(); // Any exception fails the entire regression executable.
    ++passed;
    std::cout << "PASS accepted: " << label << '\n';
}

template<class Callable>
void rejected(std::string_view label, Reason reason, std::string_view condition,
              Callable&& callable) {
    bool caught = false;
    try {
        callable();
    } catch (const checks::AcceptanceFailure& error) {
        test_assert(error.reason() == reason, "wrong acceptance failure reason");
        const std::string message = error.what();
        constexpr std::string_view prefix = "condition=";
        const auto measured = message.find("; measured[");
        test_assert(measured != std::string::npos, "missing measured delimiter");
        test_assert(message.starts_with(prefix) && measured >= prefix.size() &&
                        std::string_view(message).substr(prefix.size(), measured - prefix.size()) == condition,
                    "wrong historical condition diagnostic");
        test_assert(message.find("; measured[") != std::string::npos &&
                        message.find("; allowed[") != std::string::npos,
                    "failure diagnostic must retain measurements and limits");
        caught = true;
    }
    // Unrelated exception types propagate to main and fail, rather than count
    // as a successful rejection. A process crash also fails CTest normally.
    test_assert(caught, "invalid fixture was accepted");
    ++passed;
    std::cout << "PASS rejected: " << label << '\n';
}

void test_diagnostic_contract() {
    const auto reject_malformed = [](std::string_view label, const std::string& diagnostic,
                                     std::string_view expected_error) {
        const auto previous_passed = passed;
        bool rejected_message = false;
        try {
            rejected(label, Reason::true_residual, "expected condition", [&] {
                throw checks::AcceptanceFailure(Reason::true_residual, diagnostic);
            });
        } catch (const checks::AcceptanceFailure&) {
            throw std::runtime_error("fixture harness leaked the acceptance exception");
        } catch (const std::runtime_error& error) {
            test_assert(std::string_view(error.what()) == expected_error,
                        "malformed diagnostic failed for the wrong reason");
            rejected_message = true;
        }
        test_assert(rejected_message && passed == previous_passed,
                    "malformed diagnostic was counted as a successful acceptance fixture");
        ++passed;
        std::cout << "PASS rejected: " << label << '\n';
    };
    reject_malformed("prefix-only condition match",
                     "condition=expected condition with extra text; measured[x]=1; allowed[x]=0",
                     "wrong historical condition diagnostic");
    reject_malformed("missing measured delimiter", "condition=expected condition",
                     "missing measured delimiter");
}

// Synthetic fixtures test acceptance policy, not measured experiment output.
struct Metrics {
    GmresStatus status = GmresStatus::CONVERGED;
    long double residual = 1.0e-14L;
    long double algebraic_l2 = 1.0e-9L;
    long double algebraic_linf = 1.0e-9L;
    long double total_l2 = 0.01L;
    long double total_linf = 0.02L;
    long double observed = 2.0L;
};

void validate(const Metrics& metrics, Profile profile) {
    // These are calls to the exact policy used by both convergence drivers,
    // not a second implementation of its comparison predicates.
    checks::spectral_bounds(7.0L, 20.0L, profile);
    checks::solver_status(metrics.status, profile);
    checks::iteration_limit(2, options.max_iterations);
    checks::nonzero_rhs(1.0L, profile);
    checks::true_residual(metrics.residual, options.tolerance);
    checks::algebraic_error(metrics.algebraic_l2, metrics.algebraic_linf,
                            0.01L, 0.02L, fraction, profile);
    checks::error_budget(metrics.algebraic_l2, 1.0e-8L, 0.01L, fraction, profile);
    checks::residual_separation(metrics.residual, 0.01L, 1.0L, profile);
    checks::decreasing_error(metrics.total_l2, 0.04L, metrics.total_linf, 0.08L, profile);
    if (profile == Profile::single_mode) {
        checks::predicted_order(2.0L);
    }
    checks::observed_order(metrics.observed, 2.0L, order_slack, profile);
}

void test_profile(Profile profile) {
    const bool single = profile == Profile::single_mode;
    const std::string_view status_message = single
        ? "GMRES NOT_CONVERGED/BREAKDOWN cannot pass convergence study"
        : "NOT_CONVERGED or BREAKDOWN must fail multimode convergence";
    const std::string_view algebraic_message = single
        ? "algebraic error must be <=0.1% of discretization error in both norms"
        : "algebraic errors must be <=0.1% of spatial errors in both norms";
    const std::string_view decreasing_message = single
        ? "both spatial errors must decrease under refinement"
        : "both errors must decrease on refinement";
    const std::string_view order_message = single
        ? "observed order disagrees with analytic second-order prediction"
        : "observed order disagrees with analytic multimode prediction";
    std::cout << "PROFILE: " << (single ? "single-mode" : "multimode") << '\n';
    accepted("valid metrics", [&] { validate({}, profile); });
    for (const auto status : {GmresStatus::NOT_CONVERGED, GmresStatus::BREAKDOWN}) {
        auto bad = Metrics{};
        bad.status = status;
        rejected("nonconverged status", Reason::solver_status, status_message,
                 [&] { validate(bad, profile); });
    }
    auto bad = Metrics{};
    bad.residual = 1.0e-10L;
    rejected("true residual above tolerance", Reason::true_residual,
             "independent true residual exceeds tolerance", [&] { validate(bad, profile); });
    bad = {};
    bad.algebraic_l2 = 1.0e-4L;
    rejected("excessive algebraic L2", Reason::algebraic_error, algebraic_message,
             [&] { validate(bad, profile); });
    bad = {};
    bad.algebraic_linf = 1.0e-4L;
    rejected("excessive algebraic Linf", Reason::algebraic_error, algebraic_message,
             [&] { validate(bad, profile); });
    bad = {};
    bad.total_l2 = 0.04L;
    rejected("non-decreasing L2", Reason::non_decreasing_error, decreasing_message,
             [&] { validate(bad, profile); });
    bad = {};
    bad.total_linf = 0.08L;
    rejected("non-decreasing Linf", Reason::non_decreasing_error, decreasing_message,
             [&] { validate(bad, profile); });
    bad = {};
    bad.observed = 2.02L; // Within both broad ranges, outside the analytic slack.
    rejected("order outside analytical tolerance", Reason::observed_order, order_message,
             [&] { validate(bad, profile); });

    const auto infinity = std::numeric_limits<long double>::infinity();
    for (const auto value : {std::numeric_limits<long double>::quiet_NaN(), infinity, -infinity}) {
        rejected("nonfinite true residual", Reason::true_residual,
                 "independent true residual exceeds tolerance",
                 [&] { checks::true_residual(value, options.tolerance); });
        rejected("nonfinite observed order", Reason::observed_order, order_message,
                 [&] { checks::observed_order(value, 2.0L, order_slack, profile); });
        const std::string_view finite_message = single
            ? "nonfinite value in independent error/residual calculation"
            : "nonfinite diagnostic vector component";
        rejected("nonfinite real component", Reason::nonfinite_value, finite_message,
                 [&] { checks::finite_value({value, 0.0L}, finite_message); });
        rejected("nonfinite imaginary component", Reason::nonfinite_value, finite_message,
                 [&] { checks::finite_value({0.0L, value}, finite_message); });
    }
    for (const auto value : {std::numeric_limits<long double>::quiet_NaN(), infinity}) {
        rejected("nonfinite algebraic L2", Reason::algebraic_error, algebraic_message,
                 [&] { checks::algebraic_error(value, 0.0L, 1.0L, 1.0L, fraction, profile); });
        rejected("nonfinite algebraic Linf", Reason::algebraic_error, algebraic_message,
                 [&] { checks::algebraic_error(0.0L, value, 1.0L, 1.0L, fraction, profile); });
    }

    accepted("residual equal to tolerance", [&] {
        checks::true_residual(options.tolerance, options.tolerance);
    });
    rejected("residual nextafter tolerance", Reason::true_residual,
             "independent true residual exceeds tolerance", [&] {
        checks::true_residual(std::nextafter(static_cast<long double>(options.tolerance), infinity),
                              options.tolerance);
    });
    // Compute exactly the same inclusive threshold expression used by policy.
    const auto l2_limit = fraction * 0.01L;
    const auto inf_limit = fraction * 0.02L;
    accepted("algebraic errors equal to thresholds", [&] {
        checks::algebraic_error(l2_limit, inf_limit, 0.01L, 0.02L, fraction, profile);
    });
    rejected("algebraic L2 nextafter threshold", Reason::algebraic_error, algebraic_message, [&] {
        checks::algebraic_error(std::nextafter(l2_limit, infinity), 0.0L,
                                0.01L, 0.02L, fraction, profile);
    });
    rejected("algebraic Linf nextafter threshold", Reason::algebraic_error, algebraic_message, [&] {
        checks::algebraic_error(0.0L, std::nextafter(inf_limit, infinity),
                                0.01L, 0.02L, fraction, profile);
    });
    // A binary-exact synthetic slack isolates <= from subtraction rounding;
    // experiment drivers retain their original logarithmic slack unchanged.
    constexpr long double binary_slack = 1.0L / 1024.0L;
    const auto order_limit = 2.0L + binary_slack;
    accepted("order equal to analytical slack", [&] {
        checks::observed_order(order_limit, 2.0L, binary_slack, profile);
    });
    rejected("order nextafter analytical slack", Reason::observed_order, order_message, [&] {
        checks::observed_order(std::nextafter(order_limit, infinity), 2.0L, binary_slack, profile);
    });
    // Wide synthetic slack isolates the fixed broad interval from the separate
    // analytic-distance check; no experiment tolerance is modified.
    const long double upper = single ? 2.03L : 2.06L;
    for (const auto endpoint : {1.99L, upper}) {
        accepted("observed order at broad interval endpoint", [&] {
            checks::observed_order(endpoint, 2.0L, 0.125L, profile);
        });
    }
    for (const auto outside : {std::nextafter(1.99L, -infinity), std::nextafter(upper, infinity)}) {
        rejected("observed order outside broad interval", Reason::observed_order, order_message, [&] {
            checks::observed_order(outside, 2.0L, 0.125L, profile);
        });
    }
    if (single) {
        for (const auto endpoint : {2.0L, 2.02L}) {
            accepted("predicted order at inclusive endpoint", [&] { checks::predicted_order(endpoint); });
        }
        for (const auto outside : {std::nextafter(2.0L, -infinity),
                                    std::nextafter(2.02L, infinity),
                                    std::numeric_limits<long double>::quiet_NaN(), infinity}) {
            rejected("invalid single-mode prediction", Reason::predicted_order,
                     "analytic modal prediction must be in second-order regime",
                     [&] { checks::predicted_order(outside); });
        }
    } else {
        for (const auto outside : {std::nextafter(1.99L, -infinity),
                                    std::nextafter(upper, infinity),
                                    std::numeric_limits<long double>::quiet_NaN(), infinity}) {
            rejected("invalid multimode prediction", Reason::observed_order, order_message, [&] {
                checks::observed_order(2.0L, outside, 0.125L, profile);
            });
        }
    }
    const auto budget_message = single
        ? "residual and rounding error bound must not mask spatial error"
        : "residual/rounding-aware L2 budget must not mask spatial error";
    accepted("algebraic error and budget equal inclusive limits", [&] {
        checks::error_budget(l2_limit, l2_limit, 0.01L, fraction, profile);
    });
    rejected("budget nextafter spatial limit", Reason::error_budget, budget_message, [&] {
        checks::error_budget(0.0L, std::nextafter(l2_limit, infinity), 0.01L, fraction, profile);
    });
    const auto bounded_message = single
        ? "algebraic error exceeds residual/conditioning-aware bound"
        : "residual/rounding-aware L2 budget must not mask spatial error";
    rejected("algebraic error nextafter budget", Reason::error_budget, bounded_message, [&] {
        checks::error_budget(std::nextafter(l2_limit, infinity), l2_limit,
                             0.01L, fraction, profile);
    });
    const auto separation_message = single
        ? "residual must be at least 100 times below relative spatial error"
        : "residual must be 100 times below relative spatial error";
    const auto separation_limit = 0.01L * 0.01L / 1.0L;
    accepted("residual equal to spatial separation limit", [&] {
        checks::residual_separation(separation_limit, 0.01L, 1.0L, profile);
    });
    rejected("residual nextafter spatial separation limit", Reason::residual_separation,
             separation_message, [&] {
        checks::residual_separation(std::nextafter(separation_limit, infinity),
                                    0.01L, 1.0L, profile);
    });
    accepted("iteration count equal to global limit", [] { checks::iteration_limit(200, 200); });
    rejected("iteration count above global limit", Reason::iteration_limit,
             "GMRES global iteration limit", [] { checks::iteration_limit(201, 200); });
    rejected("zero iteration count", Reason::iteration_limit,
             "GMRES global iteration limit", [] { checks::iteration_limit(0, 200); });
    rejected("zero RHS norm", Reason::nonzero_rhs,
             single ? "nonzero continuous forcing" : "nonzero forcing",
             [&] { checks::nonzero_rhs(0.0L, profile); });
    const auto spectral_message = single
        ? "nonresonant SPD spectral bounds" : "nonresonant positive spectrum";
    accepted("spectral gap and radius at inclusive bounds", [&] {
        checks::spectral_bounds(7.0L, 7.0L, profile);
    });
    rejected("gap below spectral lower bound", Reason::spectral_bounds, spectral_message, [&] {
        checks::spectral_bounds(std::nextafter(7.0L, -infinity), 7.0L, profile);
    });
    rejected("radius below spectral gap", Reason::spectral_bounds, spectral_message, [&] {
        checks::spectral_bounds(7.0L, std::nextafter(7.0L, -infinity), profile);
    });
}

void test_limited_gmres() {
    const Grid2D grid{7, 7, 1.0, 1.0};
    const auto matrix = assemble_helmholtz_operator(grid, 1.0);
    auto rhs = continuous_manufactured_rhs(grid, 1, 1, 1.0);
    const auto second = continuous_manufactured_rhs(grid, 2, 3, 1.0, {-0.3, 0.0});
    for (std::size_t p = 0; p < rhs.size(); ++p) {
        rhs[p] += second[p];
    }
    const auto result = gmres(matrix, rhs, {options.tolerance, options.restart, 1});
    test_assert(result.status == GmresStatus::NOT_CONVERGED && result.iterations == 1,
                "real two-mode solve must reach limit without false convergence");
    const auto product = matrix.multiply(result.solution);
    long double rnorm = 0.0L;
    long double bnorm = 0.0L;
    for (std::size_t p = 0; p < rhs.size(); ++p) {
        rnorm = std::hypot(rnorm, std::abs(std::complex<long double>(rhs[p]) -
                                         std::complex<long double>(product[p])));
        bnorm = std::hypot(bnorm, std::abs(std::complex<long double>(rhs[p])));
    }
    const auto residual = rnorm / bnorm;
    test_assert(std::isfinite(residual) && residual > options.tolerance,
                "independent residual must confirm the real limited solve is inadequate");
    for (const auto profile : {Profile::single_mode, Profile::multimode}) {
        rejected("real limited-iteration result", Reason::solver_status,
                 profile == Profile::single_mode
                     ? "GMRES NOT_CONVERGED/BREAKDOWN cannot pass convergence study"
                     : "NOT_CONVERGED or BREAKDOWN must fail multimode convergence",
                 [&] { checks::solver_status(result.status, profile); });
    }
    std::cout << "limited GMRES: status=NOT_CONVERGED, iterations=" << result.iterations
              << ", independent_true_residual=" << residual << '\n';
}

} // namespace

int main() {
    try {
        test_profile(Profile::single_mode);
        test_profile(Profile::multimode);
        accepted("multiple iterations exercised", [] { checks::multiple_iterations(true); });
        rejected("multiple iterations not exercised", Reason::multiple_iterations,
                 "at least one grid must require multiple GMRES iterations",
                 [] { checks::multiple_iterations(false); });
        test_limited_gmres();
        test_diagnostic_contract();
        std::cout << "Acceptance scenarios passed: " << passed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL acceptance regression: " << error.what() << '\n';
        return 1;
    }
}
