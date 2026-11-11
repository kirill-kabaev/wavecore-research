#include "wavecore/helmholtz/assembly.hpp"
#include "wavecore/helmholtz/gmres.hpp"
#include "verification_checks.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace wavecore::helmholtz;
using Wide = std::complex<long double>;
using Reference = std::vector<Wide>;

namespace checks = wavecore::test::verification;
using checks::require;
constexpr auto profile = checks::Profile::single_mode;

struct Norms {
    long double euclidean = 0.0L;
    long double linf = 0.0L;

    void add(Wide value) {
        checks::finite_value(value, "nonfinite value in independent error/residual calculation");
        const auto magnitude = std::abs(value);
        euclidean = std::hypot(euclidean, magnitude);
        linf = std::max(linf, magnitude);
    }

    long double l2h(const Grid2D& grid) const {
        return std::sqrt(static_cast<long double>(grid.hx()) * grid.hy()) * euclidean;
    }
};

template<class Left, class Right>
Norms difference(const Left& left, const Right& right) {
    require(left.size() == right.size(), "error vector dimensions",
                {{"left.size", static_cast<long double>(left.size()), "== right.size", static_cast<long double>(right.size())}});
    Norms result;
    for (std::size_t p = 0; p < left.size(); ++p) {
        result.add(Wide(left[p]) - Wide(right[p]));
    }
    return result;
}

// Unit square in metres; TM_z, exp(-i omega t), lossless homogeneous medium,
// zero Dirichlet/PEC walls. Field amplitude is 1 V/m, k=1 inverse metre.
// Only interior nodes are stored: p=j*nx+i, x=(i+1)hx, y=(j+1)hy.
void run() {
    constexpr std::array<std::size_t, 4> sizes{7, 15, 31, 63};
    constexpr double k = 1.0;
    constexpr GmresOptions options{1.0e-11, 20, 200};
    constexpr long double pi = std::numbers::pi_v<long double>;
    constexpr long double lambda = 2.0L * pi * pi;
    constexpr long double solver_fraction = 1.0e-3L;
    constexpr long double epsilon = std::numeric_limits<double>::epsilon();

    std::cout << std::setprecision(17)
              << "format=wavecore.grid_convergence.v1, rhs=continuous, "
                 "Lx=Ly=1, k=1, m=n=1, amplitude=1, x0=0, tolerance="
              << options.tolerance << ", restart=" << options.restart
              << ", max_iterations=" << options.max_iterations << '\n'
              << "nx,ny,h,iterations,status,true_residual,L2h,Linf,p_L2h,p_Linf,"
                 "predicted_order,algebraic_L2h,algebraic_Linf,discretization_L2h,"
                 "discretization_Linf,spectral_gap,kappa2,algebraic_L2h_bound\n";

    long double previous_h = 0.0L;
    long double previous_l2 = 0.0L;
    long double previous_inf = 0.0L;
    long double previous_factor = 0.0L;
    std::size_t current_n = 0;
    long double current_h = 0.0L;
    try {
        for (const auto n : sizes) {
            current_n = n;
            const Grid2D grid{n, n, 1.0, 1.0};
            const long double h = grid.hx();
            current_h = h;
            require(grid.hx() == grid.hy(), "isotropic refinement",
                    {{"hx", grid.hx(), "== hy", grid.hy()}});
            require(previous_h == 0.0L || previous_h == 2.0L * h,
                    "each refinement must halve h exactly",
                    {{"previous_h", previous_h, "== 2*h (or zero on first grid)", 2.0L * h}});
            const auto matrix = assemble_helmholtz_operator(grid, k);
            const auto exact = manufactured_solution(grid, 1, 1);
            // Continuous forcing, NEVER a matvec of the sampled exact solution.
            const auto rhs = continuous_manufactured_rhs(grid, 1, 1, k);

            // Independent analytic eigenvalues (not the production eigenvalue helper).
            // The smallest/largest Dirichlet eigenvalues occur at (1,1)/(n,n).
            const auto sine = std::sin(pi * h / 2.0L);
            const auto cosine = std::cos(pi * h / 2.0L);
            const long double lambda_h = 8.0L * sine * sine / (h * h);
            const long double gap = lambda_h - k * k;
            const long double radius = 8.0L * cosine * cosine / (h * h) - k * k;
            // sin(z)>=2z/pi for 0<=z<=pi/2 gives lambda_min>=8 m^-2,
            // hence every discrete eigenvalue of A is >=7 m^-2 for k=1.
            checks::spectral_bounds(gap, radius, profile);
            const long double factor = (lambda - lambda_h) / gap;
            require(factor > 0.0L, "nonzero positive discretization error",
                    {{"factor", factor, ">", 0.0L}});
            Reference continuous(grid.unknown_count());
            Reference discrete(grid.unknown_count());
            Reference forcing(grid.unknown_count());
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    const auto p = grid.linear_index(i, j);
                    const long double value = std::sin(pi * grid.x(i)) *
                                              std::sin(pi * grid.y(j));
                    continuous[p] = value;
                    discrete[p] = (1.0L + factor) * value;
                    forcing[p] = (lambda - k * k) * value;
                }
            }
            require(difference(exact, continuous).linf <= 16.0L * epsilon,
                    "sampled exact mode agrees with independent coordinates",
                    {{"sample Linf difference", difference(exact, continuous).linf, "<=", 16.0L * epsilon}});
            require(difference(rhs, forcing).linf <= 64.0L * epsilon * (lambda - k * k),
                    "continuous RHS agrees with independent PDE formula",
                    {{"RHS Linf difference", difference(rhs, forcing).linf, "<=", 64.0L * epsilon * (lambda - k * k)}});

            const auto result = gmres(matrix, rhs, options);
            checks::solver_status(result.status, profile);
            checks::iteration_limit(result.iterations, options.max_iterations);
            // Fresh public CSR matvec and wide norms, independent of solver bookkeeping.
            const auto product = matrix.multiply(result.solution);
            const auto residual = difference(rhs, product);
            Norms rhs_norm;
            for (const auto value : rhs) {
                rhs_norm.add(Wide(value));
            }
            checks::nonzero_rhs(rhs_norm.euclidean, profile);
            const auto true_residual = residual.euclidean / rhs_norm.euclidean;
            checks::true_residual(true_residual, options.tolerance);

            const auto total = difference(result.solution, continuous);
            const auto algebraic = difference(result.solution, discrete);
            const auto discretization = difference(discrete, continuous);
            checks::algebraic_error(algebraic.l2h(grid), algebraic.linf,
                                    discretization.l2h(grid), discretization.linf,
                                    solver_fraction, profile);
            // ||A^-1 r||_2 <= ||r||_2/gap. Account conservatively for rounded
            // forcing, reference construction and sparse matvec using 64*eps*||A||.
            Norms reference_norm;
            Norms continuous_norm;
            for (const auto value : discrete) {
                reference_norm.add(value);
            }
            for (const auto value : continuous) {
                continuous_norm.add(value);
            }
            const long double algebraic_l2_bound = residual.l2h(grid) / gap +
                64.0L * epsilon * radius / gap * reference_norm.l2h(grid);
            checks::error_budget(algebraic.l2h(grid), algebraic_l2_bound,
                                 discretization.l2h(grid), solver_fraction, profile);

            checks::residual_separation(true_residual, discretization.euclidean,
                                        continuous_norm.euclidean, profile);
            // For odd n the centre is a node: ||u*||_inf=1, ||u*||_2,h=1/2.
            require(std::abs(discretization.linf - factor) <= 16.0L * epsilon &&
                        std::abs(discretization.l2h(grid) - factor / 2.0L) <= 16.0L * epsilon,
                    "reference norms must agree with the analytic mode error",
                    {{"Linf identity defect", std::abs(discretization.linf - factor), "<=", 16.0L * epsilon},
                     {"L2(h) identity defect", std::abs(discretization.l2h(grid) - factor / 2.0L), "<=", 16.0L * epsilon}});

            long double order_l2 = 0.0L;
            long double order_inf = 0.0L;
            long double predicted_order = 0.0L;
            if (previous_h != 0.0L) {
                checks::decreasing_error(total.l2h(grid), previous_l2, total.linf, previous_inf, profile);
                order_l2 = std::log2(previous_l2 / total.l2h(grid));
                order_inf = std::log2(previous_inf / total.linf);
                predicted_order = std::log2(previous_factor / factor);
                // Exact modal errors give p in [2,2.02] for these four grids.
                // A <=0.1% contamination at each level changes p by at most this:
                const auto order_slack = std::log2((1.0L + solver_fraction) /
                                                    (1.0L - solver_fraction));
                checks::predicted_order(predicted_order);
                for (const auto order : {order_l2, order_inf}) {
                    checks::observed_order(order, predicted_order, order_slack, profile);
                }
            }
            std::cout << n << ',' << n << ',' << h << ',' << result.iterations
                      << ",CONVERGED," << true_residual << ',' << total.l2h(grid)
                      << ',' << total.linf << ',';
            if (previous_h != 0.0L) {
                std::cout << order_l2 << ',' << order_inf << ',' << predicted_order;
            } else {
                std::cout << "NA,NA,NA";
            }
            std::cout << ',' << algebraic.l2h(grid) << ',' << algebraic.linf << ','
                      << discretization.l2h(grid) << ',' << discretization.linf << ','
                      << gap << ',' << radius / gap << ',' << algebraic_l2_bound << '\n';
            previous_h = h;
            previous_l2 = total.l2h(grid);
            previous_inf = total.linf;
            previous_factor = factor;
        }
    } catch (const std::exception& error) {
        std::ostringstream failure;
        failure << std::setprecision(17) << "n=" << current_n
                << ", h=" << current_h << " m; " << error.what();
        throw std::runtime_error(failure.str());
    }
}

} // namespace

int main() {
    try {
        run();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
