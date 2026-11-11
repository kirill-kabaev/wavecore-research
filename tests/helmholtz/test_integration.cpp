#include "wavecore/helmholtz/assembly.hpp"
#include "wavecore/helmholtz/gmres.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

using namespace wavecore::helmholtz;
using Vector = CsrMatrix::value_array;
using Wide = std::complex<long double>;

void require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

long double norm(const Vector& vector) {
    long double result = 0.0L;
    for (const auto value : vector) {
        require(std::isfinite(value.real()) && std::isfinite(value.imag()),
                "nonfinite vector component");
        result = std::hypot(result, std::abs(Wide(value)));
    }
    return result;
}

long double difference_norm(const Vector& left, const Vector& right) {
    require(left.size() == right.size(), "vector size mismatch");
    long double result = 0.0L;
    for (std::size_t i = 0; i < left.size(); ++i) {
        result = std::hypot(result, std::abs(Wide(left[i]) - Wide(right[i])));
    }
    return result;
}

struct Configuration {
    std::string_view name;
    std::size_t nx;
    std::size_t ny;
    double ly;
    double k;
    std::size_t m;
    std::size_t n;
    bool mixed;
};

// Lengths are in metres, k in inverse metres. TM_z, exp(-i omega t), PEC.
constexpr Configuration cases[]{
    {"small", 3, 2, 1.0, 3.0, 1, 2, false},
    {"medium", 15, 11, 1.0, 8.0, 1, 2, false},
    {"altered", 7, 5, 0.75, 3.0, 2, 1, false},
    {"restarted", 5, 4, 1.0, 3.0, 1, 2, true}};

void run(const Configuration& config) {
    const Grid2D grid{config.nx, config.ny, 1.0, config.ly};
    const auto matrix = assemble_helmholtz_operator(grid, config.k);
    require(matrix.rows() == grid.unknown_count() &&
                matrix.cols() == grid.unknown_count(), "matrix dimensions");

    // Independent analytical spectrum of -Delta_h on the entire Dirichlet grid.
    // This real symmetric matrix has kappa_2 = max|lambda_h-k^2|/min|...|.
    long double gap = std::numeric_limits<long double>::infinity();
    long double spectral_radius = 0.0L;
    for (std::size_t q = 1; q <= grid.ny(); ++q) {
        for (std::size_t p = 1; p <= grid.nx(); ++p) {
            const long double sx = std::sin(
                p * std::numbers::pi_v<long double> / (2.0L * (grid.nx() + 1)));
            const long double sy = std::sin(
                q * std::numbers::pi_v<long double> / (2.0L * (grid.ny() + 1)));
            const long double eigenvalue =
                4.0L * sx * sx / (static_cast<long double>(grid.hx()) * grid.hx()) +
                4.0L * sy * sy / (static_cast<long double>(grid.hy()) * grid.hy()) -
                static_cast<long double>(config.k) * config.k;
            gap = std::min(gap, std::abs(eigenvalue));
            spectral_radius = std::max(spectral_radius, std::abs(eigenvalue));
        }
    }
    const long double condition = spectral_radius / gap;
    require(std::isfinite(condition) && gap > 1.0L && condition < 1000.0L,
            "test system must be safely separated from discrete resonance");

    auto exact = manufactured_solution(grid, config.m, config.n, {1.0, 0.25});
    auto modal_rhs = algebraic_manufactured_rhs(
        grid, config.m, config.n, config.k, {1.0, 0.25});
    if (config.mixed) {
        // Three eigenvalues are needed to force GMRES(2) across a restart.
        const auto add_mode = [&](std::size_t m, std::size_t n,
                                  CsrMatrix::value_type amplitude) {
            const auto mode = manufactured_solution(grid, m, n, amplitude);
            const auto rhs = algebraic_manufactured_rhs(grid, m, n, config.k, amplitude);
            for (std::size_t i = 0; i < exact.size(); ++i) {
                exact[i] += mode[i];
                modal_rhs[i] += rhs[i];
            }
        };
        add_mode(1, 1, {0.5, -0.75});
        add_mode(3, 2, {-0.25, 0.5});
    }

    // The actual algebraic RHS is the stored CSR matvec, not continuous forcing.
    const auto rhs = matrix.multiply(exact);
    const long double rhs_norm = norm(rhs);
    const long double exact_norm = norm(exact);
    require(rhs_norm > 0.0L && exact_norm > 0.0L, "nonzero reference and RHS");
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    const auto modal_gap = difference_norm(rhs, modal_rhs) / rhs_norm;
    require(std::isfinite(modal_gap) && modal_gap <= 64.0L * epsilon * condition,
            "CSR RHS must agree with independent modal algebraic RHS");

    const GmresOptions options{1.0e-12, config.mixed ? 2U : 10U, 1000};
    const Vector initial(grid.unknown_count(), {0.0, 0.0});
    const auto result = gmres(matrix, rhs, options, initial);
    require(result.status == GmresStatus::CONVERGED,
            "GMRES must converge: BREAKDOWN and NOT_CONVERGED are failures");
    require(result.solution.size() == exact.size(), "solution dimensions");
    (void)norm(result.solution); // Reject nonfinite solutions before comparisons.

    // Independent of GMRES bookkeeping: a fresh public CSR matvec and wide norms.
    const auto product = matrix.multiply(result.solution);
    (void)norm(product);
    const long double residual = difference_norm(rhs, product) / rhs_norm;
    const long double error = difference_norm(result.solution, exact) / exact_norm;
    // Forward-error scale kappa*rho, plus a conservative O(epsilon*kappa)
    // allowance for constructing the RHS and evaluating the sparse products.
    const long double error_limit =
        2.0L * condition * (options.tolerance + 64.0L * epsilon);
    require(std::isfinite(residual) && residual <= options.tolerance,
            "independent true residual exceeds solver tolerance");
    require(std::isfinite(error) && error <= error_limit,
            "relative solution error exceeds conditioning-aware tolerance");
    require(std::isfinite(result.final_residual) &&
                std::abs(result.final_residual - residual) <= 64.0L * epsilon,
            "reported and independently computed residual disagree");
    require(result.iterations > 0 && result.iterations <= options.max_iterations,
            "global iteration count");
    require(result.residual_history.size() == result.iterations + 1,
            "residual history dimensions");
    for (const auto value : result.residual_history) {
        require(std::isfinite(value) && value >= 0.0, "finite residual history");
    }
    if (config.mixed) {
        require(result.iterations > options.restart, "restart must actually occur");
    }

    std::cout << std::setprecision(17) << config.name
              << ": nx=" << grid.nx() << ", ny=" << grid.ny()
              << ", hx=" << grid.hx() << ", hy=" << grid.hy()
              << ", k=" << config.k << ", nnz=" << matrix.nnz()
              << ", spectral_gap=" << gap << ", kappa2=" << condition
              << ", restart=" << options.restart
              << ", max_iterations=" << options.max_iterations
              << ", tolerance=" << options.tolerance
              << ", status=CONVERGED, iterations=" << result.iterations
              << ", true_relative_residual=" << residual
              << ", relative_solution_error=" << error
              << ", error_limit=" << error_limit << '\n';
    std::cout << "residual_history:";
    for (const auto value : result.residual_history) {
        std::cout << ' ' << value;
    }
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: wavecore_integration_test CASE");
        for (const auto& config : cases) {
            if (config.name == argv[1]) {
                run(config);
                return 0;
            }
        }
        throw std::runtime_error("unknown integration case");
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
