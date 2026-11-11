#include "wavecore/helmholtz/assembly.hpp"
#include "wavecore/helmholtz/gmres.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <initializer_list>
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

struct Measurement {
    std::string_view name;
    long double value;
    std::string_view relation;
    long double limit;
};

void require(bool condition, std::string_view message,
             std::initializer_list<Measurement> measurements = {}) {
    if (!condition) {
        std::ostringstream failure;
        failure << std::setprecision(17) << "condition=" << message;
        for (const auto& item : measurements) {
            failure << "; measured[" << item.name << "]=" << item.value
                    << "; allowed[" << item.name << "]=" << item.relation
                    << ' ' << item.limit;
        }
        throw std::runtime_error(failure.str());
    }
}

struct Norms {
    long double l2 = 0.0L;
    long double linf = 0.0L;

    void add(Wide value) {
        require(std::isfinite(value.real()) && std::isfinite(value.imag()),
                "nonfinite diagnostic vector component",
                {{"abs(real)", std::abs(value.real()), "<= finite maximum", std::numeric_limits<long double>::max()},
                 {"abs(imag)", std::abs(value.imag()), "<= finite maximum", std::numeric_limits<long double>::max()}});
        const auto magnitude = std::abs(value);
        l2 = std::hypot(l2, magnitude);
        linf = std::max(linf, magnitude);
    }

    long double l2h(long double h) const { return h * l2; }
};

template<class Left, class Right>
Norms difference(const Left& left, const Right& right) {
    require(left.size() == right.size(), "diagnostic vector dimensions",
                {{"left.size", static_cast<long double>(left.size()), "== right.size", static_cast<long double>(right.size())}});
    Norms result;
    for (std::size_t p = 0; p < left.size(); ++p) {
        result.add(Wide(left[p]) - Wide(right[p]));
    }
    return result;
}

// Independent references for diagnostics only; production assembly and
// manufactured data are always obtained through the existing public API.
long double modal_eigenvalue(long double h, unsigned m, unsigned n) {
    constexpr auto pi = std::numbers::pi_v<long double>;
    const auto sx = std::sin(m * pi * h / 2.0L);
    const auto sy = std::sin(n * pi * h / 2.0L);
    return 4.0L * (sx * sx + sy * sy) / (h * h);
}

// TM_z, exp(-i omega t), homogeneous lossless medium, PEC/zero Dirichlet,
// unit square in metres, k=1 m^-1. Amplitudes 1 and -0.3 V/m.
// Interior nodes only, nx=ny=n, h=1/(n+1), p=j*n+i (x fastest).
void run() {
    constexpr std::array<std::size_t, 4> sizes{7, 15, 31, 63};
    constexpr GmresOptions options{1.0e-11, 20, 200};
    constexpr long double pi = std::numbers::pi_v<long double>;
    constexpr long double lambda1 = 2.0L * pi * pi;
    constexpr long double lambda2 = 13.0L * pi * pi;
    constexpr long double weight = -0.3L;
    constexpr long double epsilon = std::numeric_limits<double>::epsilon();
    constexpr long double solver_fraction = 0.001L;
    const auto order_slack = std::log2((1.0L + solver_fraction) /
                                     (1.0L - solver_fraction));

    std::cout << std::setprecision(17)
              << "format=wavecore.multimode_convergence.v1, rhs=continuous, "
                 "Lx=Ly=1, k=1, modes=(1;1)+(2;3), amplitudes=1;-0.3, "
                 "x0=0, tolerance=" << options.tolerance
              << ", restart=" << options.restart
              << ", max_iterations=" << options.max_iterations << '\n'
              << "nx,ny,h,iterations,status,true_residual,L2h,Linf,p_L2h,p_Linf,"
                 "predicted_p_L2h,predicted_p_Linf,algebraic_L2h,algebraic_Linf,"
                 "discretization_L2h,discretization_Linf,spectral_gap,kappa2,"
                 "algebraic_L2h_bound\n";

    long double previous_h = 0.0L;
    Norms previous_total;
    Norms previous_spatial;
    bool multiple_iterations = false;
    std::size_t current_n = 0;
    long double current_h = 0.0L;
    try {
        for (const auto n : sizes) {
            current_n = n;
            const Grid2D grid{n, n, 1.0, 1.0};
            const long double h = grid.hx();
            current_h = h;
            require(grid.hx() == grid.hy() &&
                        (previous_h == 0.0L || previous_h == 2.0L * h),
                    "isotropic grid spacing must halve at each level",
                    {{"hx", grid.hx(), "== hy", grid.hy()},
                     {"previous_h", previous_h, "== 2*h (or zero on first grid)", 2.0L * h}});
            const auto matrix = assemble_helmholtz_operator(grid, 1.0);
            auto sampled = manufactured_solution(grid, 1, 1);
            const auto second = manufactured_solution(grid, 2, 3, {-0.3, 0.0});
            auto rhs = continuous_manufactured_rhs(grid, 1, 1, 1.0);
            const auto rhs_second = continuous_manufactured_rhs(grid, 2, 3, 1.0,
                                                                {-0.3, 0.0});
            for (std::size_t p = 0; p < rhs.size(); ++p) {
                sampled[p] += second[p];
                rhs[p] += rhs_second[p]; // Continuous forcing, never A_h*u_exact.
            }

            const auto gap = modal_eigenvalue(h, 1, 1) - 1.0L;
            const auto radius = modal_eigenvalue(h, static_cast<unsigned>(n),
                                                static_cast<unsigned>(n)) - 1.0L;
            // sin(z)>=2z/pi implies lambda_min(-Delta_h)>=8, hence gap>=7.
            require(gap >= 7.0L && radius >= gap, "nonresonant positive spectrum",
                    {{"gap", gap, ">=", 7.0L}, {"radius", radius, ">= gap", gap}});
            const auto d1 = (lambda1 - modal_eigenvalue(h, 1, 1)) / gap;
            const auto d2 = (lambda2 - modal_eigenvalue(h, 2, 3)) /
                            (modal_eigenvalue(h, 2, 3) - 1.0L);
            require(d1 > 0.0L && d2 > 0.0L, "positive modal spatial error factors",
                    {{"d1", d1, ">", 0.0L}, {"d2", d2, ">", 0.0L}});
            Reference continuous(grid.unknown_count());
            Reference discrete(grid.unknown_count());
            Reference forcing(grid.unknown_count());
            Norms continuous_norm;
            Norms discrete_norm;
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    const auto p = grid.linear_index(i, j);
                    const long double x = grid.x(i);
                    const long double y = grid.y(j);
                    const auto phi1 = std::sin(pi * x) * std::sin(pi * y);
                    const auto phi2 = std::sin(2.0L * pi * x) * std::sin(3.0L * pi * y);
                    continuous[p] = phi1 + weight * phi2;
                    discrete[p] = (1.0L + d1) * phi1 + weight * (1.0L + d2) * phi2;
                    forcing[p] = (lambda1 - 1.0L) * phi1 + weight * (lambda2 - 1.0L) * phi2;
                    continuous_norm.add(continuous[p]);
                    discrete_norm.add(discrete[p]);
                }
            }
            require(difference(sampled, continuous).linf <= 64.0L * epsilon,
                    "manufactured mixture must match independent coordinates/formula",
                    {{"sample Linf difference", difference(sampled, continuous).linf, "<=", 64.0L * epsilon}});
            require(difference(rhs, forcing).linf <= 64.0L * epsilon *
                        (lambda1 - 1.0L + std::abs(weight) * (lambda2 - 1.0L)),
                    "RHS must match the independent continuous PDE formula",
                    {{"RHS Linf difference", difference(rhs, forcing).linf, "<=",
                      64.0L * epsilon * (lambda1 - 1.0L + std::abs(weight) * (lambda2 - 1.0L))}});

            const auto result = gmres(matrix, rhs, options);
            require(result.status == GmresStatus::CONVERGED,
                    "NOT_CONVERGED or BREAKDOWN must fail multimode convergence",
                    {{"status (0=CONVERGED,1=NOT_CONVERGED,2=BREAKDOWN)", static_cast<long double>(result.status), "== CONVERGED", static_cast<long double>(GmresStatus::CONVERGED)}});
            require(result.iterations > 0 && result.iterations <= options.max_iterations,
                    "GMRES global iteration limit",
                    {{"iterations", static_cast<long double>(result.iterations), ">", 0.0L},
                     {"iterations", static_cast<long double>(result.iterations), "<= max_iterations", static_cast<long double>(options.max_iterations)}});
            multiple_iterations = multiple_iterations || result.iterations > 1;
            // Fresh public CSR product; no use of reported residual for acceptance.
            const auto product = matrix.multiply(result.solution);
            const auto residual = difference(rhs, product);
            Norms rhs_norm;
            for (const auto value : rhs) {
                rhs_norm.add(Wide(value));
            }
            require(rhs_norm.l2 > 0.0L, "nonzero forcing",
                    {{"RHS norm", rhs_norm.l2, ">", 0.0L}});
            const auto true_residual = residual.l2 / rhs_norm.l2;
            require(std::isfinite(true_residual) && true_residual <= options.tolerance,
                    "independent true residual exceeds tolerance",
                    {{"true_residual", true_residual, "<= tolerance", options.tolerance}});

            const auto total = difference(result.solution, continuous);
            const auto algebraic = difference(result.solution, discrete);
            const auto spatial = difference(discrete, continuous);
            // Discrete sine orthogonality gives an independent closed-form L2 error.
            const auto analytic_l2 = 0.5L * std::hypot(d1, weight * d2);
            require(std::abs(spatial.l2h(h) - analytic_l2) <= 16.0L * epsilon,
                    "spatial L2 reference must agree with orthogonal mode formula",
                    {{"L2(h) identity defect", std::abs(spatial.l2h(h) - analytic_l2), "<=", 16.0L * epsilon}});
            require(algebraic.l2 <= solver_fraction * spatial.l2 &&
                        algebraic.linf <= solver_fraction * spatial.linf,
                    "algebraic errors must be <=0.1% of spatial errors in both norms",
                    {{"algebraic L2", algebraic.l2, "<=", solver_fraction * spatial.l2},
                     {"algebraic Linf", algebraic.linf, "<=", solver_fraction * spatial.linf}});
            const auto algebraic_bound = residual.l2h(h) / gap +
                64.0L * epsilon * (radius / gap) * discrete_norm.l2h(h);
            require(algebraic.l2h(h) <= algebraic_bound &&
                        algebraic_bound <= solver_fraction * spatial.l2h(h),
                    "residual/rounding-aware L2 budget must not mask spatial error",
                    {{"algebraic L2(h)", algebraic.l2h(h), "<= error budget", algebraic_bound},
                     {"L2 error budget", algebraic_bound, "<=", solver_fraction * spatial.l2h(h)}});
            require(true_residual <= 0.01L * spatial.l2 / continuous_norm.l2,
                    "residual must be 100 times below relative spatial error",
                    {{"true_residual", true_residual, "<=", 0.01L * spatial.l2 / continuous_norm.l2}});

            std::cout << n << ',' << n << ',' << h << ',' << result.iterations
                      << ",CONVERGED," << true_residual << ',' << total.l2h(h)
                      << ',' << total.linf << ',';
            if (previous_h == 0.0L) {
                std::cout << "NA,NA,NA,NA";
            } else {
                require(total.l2h(h) < previous_total.l2h(previous_h) &&
                            total.linf < previous_total.linf,
                        "both errors must decrease on refinement",
                    {{"L2(h)", total.l2h(h), "< previous L2(h)", previous_total.l2h(previous_h)},
                         {"Linf", total.linf, "< previous Linf", previous_total.linf}});
                const auto p_l2 = std::log2(previous_total.l2h(previous_h) / total.l2h(h));
                const auto p_inf = std::log2(previous_total.linf / total.linf);
                const auto predicted_l2 = std::log2(previous_spatial.l2h(previous_h) /
                                                    spatial.l2h(h));
                const auto predicted_inf = std::log2(previous_spatial.linf / spatial.linf);
                // Analytic mixture predictions are approx. 2.0514/2.0477,
                // 2.0126/1.9993, 2.0031/2.0029 (L2/Linf). Nodal maxima can shift.
                // A <=0.1% solver contribution perturbs p by <=order_slack.
                const auto check_order = [&](long double observed, long double predicted) {
                    require(std::isfinite(observed) && predicted >= 1.99L &&
                                predicted <= 2.06L && observed >= 1.99L && observed <= 2.06L &&
                                std::abs(observed - predicted) <= order_slack,
                            "observed order disagrees with analytic multimode prediction",
                    {{"finite(observed)", std::isfinite(observed) ? 1.0L : 0.0L, "==", 1.0L},
                             {"predicted", predicted, ">=", 1.99L}, {"predicted", predicted, "<=", 2.06L},
                             {"observed", observed, ">=", 1.99L}, {"observed", observed, "<=", 2.06L},
                             {"abs(observed-predicted)", std::abs(observed - predicted), "<=", order_slack}});
                };
                check_order(p_l2, predicted_l2);
                check_order(p_inf, predicted_inf);
                std::cout << p_l2 << ',' << p_inf << ',' << predicted_l2 << ',' << predicted_inf;
            }
            std::cout << ',' << algebraic.l2h(h) << ',' << algebraic.linf << ','
                      << spatial.l2h(h) << ',' << spatial.linf << ',' << gap << ','
                      << radius / gap << ',' << algebraic_bound << '\n';
            previous_h = h;
            previous_total = total;
            previous_spatial = spatial;
        }
        require(multiple_iterations, "at least one grid must require multiple GMRES iterations",
                    {{"multiple_iterations", multiple_iterations ? 1.0L : 0.0L, "==", 1.0L}});
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
