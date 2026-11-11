#pragma once

#include "wavecore/helmholtz/gmres.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

// Test-only acceptance policy. Production APIs, references, norms, and solver
// arithmetic are deliberately outside this header. Keep historical diagnostics.
//
// Preconditions for normal experiment calls (not deliberately invalid fixtures):
// - Residuals, error norms, reference norms, budgets, spectral gap/radius,
//   observed/predicted orders, tolerances, fractions, and slack must be finite.
// - Norms, residuals, budgets, fraction, and slack must be nonnegative. Norms
//   used as denominators must be strictly positive; Stage 5B also has strictly
//   positive spatial errors. Tolerance and the global iteration limit are positive.
// - finite_value validates real/imaginary components BEFORE norm accumulation.
//   Callers must also use representable norm arithmetic: finite components alone
//   do not guarantee a finite aggregate norm for arbitrary data.
// - nonzero_rhs assumes a previously validated RHS norm; it tests positivity,
//   not finiteness. algebraic_error, error_budget, residual_separation, and
//   decreasing_error assume validated finite, nonnegative norm inputs.
// - true_residual explicitly rejects a nonfinite residual but assumes a valid
//   tolerance. Order checks assume valid finite slack; single-mode prediction
//   validation is a separate prerequisite for the observed-order check.
//
// Driver order is significant for the first reported failure: validate finite
// samples/norm inputs and the case spectrum, then solver status and iteration
// count; validate a nonzero RHS before normalization; check true residual before
// algebraic errors, budget, and residual/spatial-error separation. On refinement,
// check decreasing errors before order. Single-mode checks predicted_order before
// observed_order; multimode checks each prediction inside observed_order.
// error_budget preserves two ordered single-mode failures and one combined
// multimode failure. These are the existing Stage 5B predicates, not general
// input validators: documenting preconditions must not add guards or change
// thresholds. The spectral bound is specific to the unit-square, k=1 cases.
namespace wavecore::test::verification {

enum class Profile { single_mode, multimode };
enum class Reason {
    other, nonfinite_value, spectral_bounds, solver_status, iteration_limit, nonzero_rhs, true_residual,
    algebraic_error, error_budget, residual_separation, non_decreasing_error,
    predicted_order, observed_order, multiple_iterations
};

class AcceptanceFailure : public std::runtime_error {
public:
    AcceptanceFailure(Reason reason, const std::string& message)
        : std::runtime_error(message), reason_(reason) {}
    [[nodiscard]] Reason reason() const noexcept { return reason_; }
private:
    Reason reason_;
};

struct Measurement {
    std::string_view name;
    long double value;
    std::string_view relation;
    long double limit;
};

inline void require(bool condition, std::string_view message,
                    std::initializer_list<Measurement> measurements = {},
                    Reason reason = Reason::other) {
    if (!condition) {
        std::ostringstream failure;
        failure << std::setprecision(17) << "condition=" << message;
        for (const auto& item : measurements) {
            failure << "; measured[" << item.name << "]=" << item.value
                    << "; allowed[" << item.name << "]=" << item.relation
                    << ' ' << item.limit;
        }
        throw AcceptanceFailure(reason, failure.str());
    }
}

inline void finite_value(std::complex<long double> value, std::string_view message) {
    require(std::isfinite(value.real()) && std::isfinite(value.imag()), message,
            {{"abs(real)", std::abs(value.real()), "<= finite maximum", std::numeric_limits<long double>::max()},
             {"abs(imag)", std::abs(value.imag()), "<= finite maximum", std::numeric_limits<long double>::max()}},
            Reason::nonfinite_value);
}

inline void spectral_bounds(long double gap, long double radius, Profile profile) {
    require(gap >= 7.0L && radius >= gap, profile == Profile::single_mode
                ? "nonresonant SPD spectral bounds" : "nonresonant positive spectrum",
            {{"gap", gap, ">=", 7.0L}, {"radius", radius, ">= gap", gap}}, Reason::spectral_bounds);
}

inline void solver_status(helmholtz::GmresStatus status, Profile profile) {
    require(status == helmholtz::GmresStatus::CONVERGED,
            profile == Profile::single_mode
                ? "GMRES NOT_CONVERGED/BREAKDOWN cannot pass convergence study"
                : "NOT_CONVERGED or BREAKDOWN must fail multimode convergence",
            {{"status (0=CONVERGED,1=NOT_CONVERGED,2=BREAKDOWN)", static_cast<long double>(status), "== CONVERGED", static_cast<long double>(helmholtz::GmresStatus::CONVERGED)}},
            Reason::solver_status);
}

inline void iteration_limit(std::size_t iterations, std::size_t max_iterations) {
    require(iterations > 0 && iterations <= max_iterations, "GMRES global iteration limit",
            {{"iterations", static_cast<long double>(iterations), ">", 0.0L},
             {"iterations", static_cast<long double>(iterations), "<= max_iterations", static_cast<long double>(max_iterations)}},
            Reason::iteration_limit);
}

// This argument must be independently computed by the driver, not copied from
// GmresResult::final_residual or an Arnoldi/QR estimate.
inline void true_residual(long double residual, long double tolerance) {
    require(std::isfinite(residual) && residual <= tolerance,
            "independent true residual exceeds tolerance",
            {{"true_residual", residual, "<= tolerance", tolerance}}, Reason::true_residual);
}

inline void nonzero_rhs(long double norm, Profile profile) {
    require(norm > 0.0L, profile == Profile::single_mode
                ? "nonzero continuous forcing" : "nonzero forcing",
            {{"RHS norm", norm, ">", 0.0L}}, Reason::nonzero_rhs);
}

// Either weighted or unweighted L2 is valid here if both errors use the same
// convention; preserve each historical driver's arithmetic and labels.
inline void algebraic_error(long double l2, long double linf,
                            long double spatial_l2, long double spatial_linf,
                            long double fraction, Profile profile) {
    require(l2 <= fraction * spatial_l2 && linf <= fraction * spatial_linf,
            profile == Profile::single_mode
                ? "algebraic error must be <=0.1% of discretization error in both norms"
                : "algebraic errors must be <=0.1% of spatial errors in both norms",
            {{profile == Profile::single_mode ? "algebraic L2(h)" : "algebraic L2", l2, "<=", fraction * spatial_l2},
             {"algebraic Linf", linf, "<=", fraction * spatial_linf}}, Reason::algebraic_error);
}

inline void error_budget(long double l2h, long double budget,
                         long double spatial_l2h, long double fraction, Profile profile) {
    const bool bounded = l2h <= budget;
    const bool separated = budget <= fraction * spatial_l2h;
    if (profile == Profile::single_mode) {
        require(bounded, "algebraic error exceeds residual/conditioning-aware bound",
                {{"algebraic L2(h)", l2h, "<= error budget", budget}}, Reason::error_budget);
        require(separated, "residual and rounding error bound must not mask spatial error",
                {{"L2 error budget", budget, "<=", fraction * spatial_l2h}}, Reason::error_budget);
    } else {
        require(bounded && separated, "residual/rounding-aware L2 budget must not mask spatial error",
                {{"algebraic L2(h)", l2h, "<= error budget", budget},
                 {"L2 error budget", budget, "<=", fraction * spatial_l2h}}, Reason::error_budget);
    }
}

inline void residual_separation(long double residual, long double spatial_l2,
                                long double continuous_l2, Profile profile) {
    require(residual <= 0.01L * spatial_l2 / continuous_l2,
            profile == Profile::single_mode
                ? "residual must be at least 100 times below relative spatial error"
                : "residual must be 100 times below relative spatial error",
            {{"true_residual", residual, "<=", 0.01L * spatial_l2 / continuous_l2}},
            Reason::residual_separation);
}

inline void decreasing_error(long double l2h, long double previous_l2h,
                              long double linf, long double previous_linf, Profile profile) {
    require(l2h < previous_l2h && linf < previous_linf,
            profile == Profile::single_mode
                ? "both spatial errors must decrease under refinement"
                : "both errors must decrease on refinement",
            {{"L2(h)", l2h, "< previous L2(h)", previous_l2h},
             {"Linf", linf, "< previous Linf", previous_linf}}, Reason::non_decreasing_error);
}

// Single-mode historically checks the prediction separately; multimode checks
// each prediction together with the corresponding observed order.
inline void predicted_order(long double predicted) {
    require(predicted >= 2.0L && predicted <= 2.02L,
            "analytic modal prediction must be in second-order regime",
            {{"predicted order", predicted, ">=", 2.0L},
             {"predicted order", predicted, "<=", 2.02L}}, Reason::predicted_order);
}

inline void observed_order(long double observed, long double predicted,
                            long double slack, Profile profile) {
    const long double upper = profile == Profile::single_mode ? 2.03L : 2.06L;
    const bool acceptable = std::isfinite(observed) &&
        (profile == Profile::single_mode || (predicted >= 1.99L && predicted <= upper)) &&
        observed >= 1.99L && observed <= upper && std::abs(observed - predicted) <= slack;
    if (profile == Profile::single_mode) {
        require(acceptable, "observed order disagrees with analytic second-order prediction",
                {{"finite(order)", std::isfinite(observed) ? 1.0L : 0.0L, "==", 1.0L},
                 {"order", observed, ">=", 1.99L}, {"order", observed, "<=", 2.03L},
                 {"abs(order-predicted)", std::abs(observed - predicted), "<=", slack}},
                Reason::observed_order);
    } else {
        require(acceptable, "observed order disagrees with analytic multimode prediction",
                {{"finite(observed)", std::isfinite(observed) ? 1.0L : 0.0L, "==", 1.0L},
                 {"predicted", predicted, ">=", 1.99L}, {"predicted", predicted, "<=", 2.06L},
                 {"observed", observed, ">=", 1.99L}, {"observed", observed, "<=", 2.06L},
                 {"abs(observed-predicted)", std::abs(observed - predicted), "<=", slack}},
                Reason::observed_order);
    }
}

inline void multiple_iterations(bool exercised) {
    require(exercised, "at least one grid must require multiple GMRES iterations",
            {{"multiple_iterations", exercised ? 1.0L : 0.0L, "==", 1.0L}},
            Reason::multiple_iterations);
}

} // namespace wavecore::test::verification
