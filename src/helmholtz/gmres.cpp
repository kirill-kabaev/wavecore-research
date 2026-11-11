#include "wavecore/helmholtz/gmres.hpp"
#include "gmres_diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace wavecore::helmholtz {
namespace {

using Complex = CsrMatrix::value_type;
using Vector = CsrMatrix::value_array;

bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void require_finite(const Vector& values) {
    for (const auto value : values) {
        if (!finite(value)) {
            throw std::overflow_error("GMRES arithmetic produced a nonfinite value");
        }
    }
}

double norm(const Vector& values) {
    double result = 0.0;
    for (const auto value : values) {
        result = std::hypot(result, value.real());
        result = std::hypot(result, value.imag());
    }
    if (!std::isfinite(result)) {
        throw std::overflow_error("GMRES vector norm is not finite");
    }
    return result;
}

Complex inner_product(const Vector& left, const Vector& right) {
    Complex result{0.0, 0.0};
    for (std::size_t i = 0; i < left.size(); ++i) {
        result += std::conj(left[i]) * right[i];
    }
    if (!finite(result)) {
        throw std::overflow_error("GMRES inner product is not finite");
    }
    return result;
}

Vector residual(const CsrMatrix& matrix, const Vector& rhs, const Vector& solution) {
    auto result = matrix.multiply(solution);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = rhs[i] - result[i];
    }
    require_finite(result);
    return result;
}

struct Rotation {
    double cosine;
    Complex sine;
};

// G = [c s; -conj(s) c], with real c, annihilates the second entry.
Rotation make_rotation(Complex a, Complex b) {
    const double abs_a = std::abs(a);
    const double abs_b = std::abs(b);
    const double magnitude = std::hypot(abs_a, abs_b);
    if (!std::isfinite(magnitude)) {
        throw std::overflow_error("GMRES Givens rotation magnitude is not finite");
    }
    if (abs_b == 0.0) {
        return {1.0, {0.0, 0.0}};
    }
    if (abs_a == 0.0) {
        return {0.0, std::conj(b) / abs_b};
    }
    return {abs_a / magnitude, (a / abs_a) * (std::conj(b) / magnitude)};
}

void rotate(const Rotation& rotation, Complex& first, Complex& second) {
    const auto old_first = first;
    first = rotation.cosine * old_first + rotation.sine * second;
    second = -std::conj(rotation.sine) * old_first + rotation.cosine * second;
    if (!finite(first) || !finite(second)) {
        throw std::overflow_error("GMRES rotated entry is not finite");
    }
}

}  // namespace

GmresResult detail::gmres_observed(const CsrMatrix& matrix,
                  const Vector& rhs,
                  const GmresOptions& options,
                  const Vector& initial_guess,
                  const ArnoldiObserver& observer) {
    const auto n = matrix.rows();
    if (matrix.cols() != n || rhs.size() != n ||
        (!initial_guess.empty() && initial_guess.size() != n)) {
        throw std::invalid_argument("GMRES requires square A and matching b and x0 sizes");
    }
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0.0 ||
        options.restart == 0 || options.max_iterations == 0) {
        throw std::invalid_argument("GMRES tolerance, restart and iteration limit must be positive");
    }
    for (const auto* vector : {&rhs, &initial_guess}) {
        for (const auto value : *vector) {
            if (!finite(value)) {
                throw std::invalid_argument("GMRES inputs must be finite");
            }
        }
    }

    GmresResult result;
    result.solution = initial_guess.empty() ? Vector(n, Complex{}) : initial_guess;
    const double rhs_norm = norm(rhs);
    const double denominator = rhs_norm == 0.0 ? 1.0 : rhs_norm;
    auto r = residual(matrix, rhs, result.solution);
    auto record_residual = [&] {
        result.final_residual = norm(r) / denominator;
        if (!std::isfinite(result.final_residual)) {
            throw std::overflow_error("GMRES normalized residual is not finite");
        }
        result.residual_history.push_back(result.final_residual);
    };
    record_residual();

    if (result.final_residual <= options.tolerance) {
        result.status = GmresStatus::CONVERGED;
        return result;
    }

    while (result.iterations < options.max_iterations) {
        const auto cycle_size = std::min({options.restart, n,
                                         options.max_iterations - result.iterations});
        if (cycle_size == std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("GMRES Krylov workspace size overflow");
        }
        const auto base = result.solution;
        const double beta = norm(r);
        std::vector<Vector> basis;
        Vector first = r;
        for (auto& value : first) {
            value /= beta;
        }
        basis.push_back(std::move(first));
        // Hessenberg columns are stored individually: h[column][row].
        std::vector<Vector> h;
        std::vector<Rotation> rotations;
        Vector g(cycle_size + 1, Complex{});
        g[0] = beta;

        for (std::size_t column = 0; column < cycle_size; ++column) {
            auto w = matrix.multiply(basis[column]);
            require_finite(w);
            const double original_norm = norm(w);
            h.emplace_back(column + 2, Complex{});
            auto& current = h.back();

            // Two MGS passes, each using <v,w> = sum conj(v_i)*w_i.
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t row = 0; row <= column; ++row) {
                    const auto projection = inner_product(basis[row], w);
                    current[row] += projection;
                    for (std::size_t i = 0; i < n; ++i) {
                        w[i] -= projection * basis[row][i];
                    }
                }
            }
            require_finite(w);
            const double next_norm = norm(w);
            const bool breakdown = next_norm <=
                16.0 * std::numeric_limits<double>::epsilon() * original_norm;
            current[column + 1] = next_norm;
            for (std::size_t previous = 0; previous < column; ++previous) {
                rotate(rotations[previous], current[previous], current[previous + 1]);
            }
            rotations.push_back(make_rotation(current[column], current[column + 1]));
            rotate(rotations.back(), current[column], current[column + 1]);
            rotate(rotations.back(), g[column], g[column + 1]);
            ++result.iterations;

            Vector coefficients(column + 1, Complex{});
            for (std::size_t reverse = column + 1; reverse > 0; --reverse) {
                const auto row = reverse - 1;
                if (h[row][row] == Complex{}) {
                    record_residual();
                    result.status = GmresStatus::BREAKDOWN;
                    return result;
                }
                Complex value = g[row];
                for (std::size_t k = row + 1; k <= column; ++k) {
                    value -= h[k][row] * coefficients[k];
                }
                coefficients[row] = value / h[row][row];
            }
            require_finite(coefficients);
            result.solution = base;
            for (std::size_t k = 0; k <= column; ++k) {
                for (std::size_t i = 0; i < n; ++i) {
                    result.solution[i] += basis[k][i] * coefficients[k];
                }
            }
            require_finite(result.solution);
            r = residual(matrix, rhs, result.solution);
            record_residual();
            if (observer) {
                observer(basis, w, next_norm, breakdown,
                         std::abs(g[column + 1]) / denominator, result.final_residual);
            }
            if (result.final_residual <= options.tolerance) {
                result.status = GmresStatus::CONVERGED;
                return result;
            }
            if (breakdown) {
                result.status = GmresStatus::BREAKDOWN;
                return result;
            }
            if (column + 1 < cycle_size) {
                for (auto& value : w) {
                    value /= next_norm;
                }
                basis.push_back(std::move(w));
            }
        }
    }
    return result;
}

GmresResult gmres(const CsrMatrix& matrix,
                  const Vector& rhs,
                  const GmresOptions& options,
                  const Vector& initial_guess) {
    return detail::gmres_observed(matrix, rhs, options, initial_guess, {});
}

}  // namespace wavecore::helmholtz
