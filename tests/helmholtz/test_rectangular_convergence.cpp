#include "verification_checks.hpp"
#include "wavecore/helmholtz/assembly.hpp"

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
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace wavecore::helmholtz;
namespace checks = wavecore::test::verification;
using checks::require;
using Wide = std::complex<long double>;
using Reference = std::vector<Wide>;
constexpr auto profile = checks::Profile::single_mode;

struct Norms {
    long double l2 = 0.0L;
    long double linf = 0.0L;
    void add(Wide value) {
        checks::finite_value(value, "finite rectangular diagnostic component");
        const auto magnitude = std::abs(value);
        l2 = std::hypot(l2, magnitude);
        linf = std::max(linf, magnitude);
    }
    long double weighted(long double weight) const { return weight * l2; }
};

template<class Left, class Right>
Norms difference(const Left& left, const Right& right) {
    require(left.size() == right.size(), "rectangular diagnostic dimensions",
            {{"left.size", static_cast<long double>(left.size()), "== right.size",
              static_cast<long double>(right.size())}});
    Norms result;
    for (std::size_t p = 0; p < left.size(); ++p) {
        result.add(Wide(left[p]) - Wide(right[p]));
    }
    return result;
}

// TM_z, exp(-i omega t), homogeneous isotropic lossless medium, PEC/zero
// Dirichlet. Lengths in metres, k=1 m^-1, field amplitude=1 V/m.
// nx/ny count interior nodes; storage p=j*nx+i is x-fastest.
void run() {
    constexpr std::array<std::array<std::size_t, 2>, 4> grids{{
        {7, 5}, {15, 11}, {31, 23}, {63, 47}}};
    constexpr long double lx = 1.0L;
    constexpr long double ly = 1.5L;
    constexpr long double k = 1.0L;
    constexpr long double pi = std::numbers::pi_v<long double>;
    constexpr long double lambda = (pi / lx) * (pi / lx) + (pi / ly) * (pi / ly);
    constexpr long double gap_lower_bound = 4.0L / (lx * lx) + 4.0L / (ly * ly) - k * k;
    constexpr long double epsilon = std::numeric_limits<double>::epsilon();
    constexpr long double fraction = 0.001L;
    constexpr GmresOptions options{1.0e-11, 20, 200};
    const auto order_slack = std::log2((1.0L + fraction) / (1.0L - fraction));
    const auto exact_continuous_l2 = std::sqrt(lx * ly) / 2.0L;

    // Pre-run analytic predictions: d_h=(lambda-lambda_h)/(lambda_h-1).
    // Orders are 2.011105127551421, 2.002768065401026, 2.000691504867757.
    // Thus the unchanged single-mode order policy [2,2.02]/[1.99,2.03]
    // and log2(1.001/0.999) contamination allowance apply to this rectangle.
    std::cout << std::setprecision(17)
              << "format=wavecore.rectangular_convergence.v1, rhs=continuous, "
                 "Lx=1, Ly=1.5, k=1, m=n=1, amplitude=1, x0=0, tolerance="
              << options.tolerance << ", restart=" << options.restart
              << ", max_iterations=" << options.max_iterations << '\n'
              << "nx,ny,hx,hy,h,iterations,status,true_residual,csr_true_residual,"
                 "L2h,Linf,p_L2h,p_Linf,predicted_order,algebraic_L2h,algebraic_Linf,"
                 "spatial_L2h,spatial_Linf,lambda,lambda_h,spectral_gap,kappa2,algebraic_L2h_bound\n";

    long double previous_hx = 0.0L;
    long double previous_hy = 0.0L;
    long double previous_factor = 0.0L;
    long double previous_total_l2 = 0.0L;
    long double previous_total_inf = 0.0L;
    std::size_t current_nx = 0;
    std::size_t current_ny = 0;
    long double current_hx = 0.0L;
    long double current_hy = 0.0L;
    try {
        for (const auto& dimensions : grids) {
            const auto nx = dimensions[0];
            const auto ny = dimensions[1];
            current_nx = nx;
            current_ny = ny;
            const Grid2D grid{nx, ny, static_cast<double>(lx), static_cast<double>(ly)};
            const long double hx = grid.hx();
            const long double hy = grid.hy();
            current_hx = hx;
            current_hy = hy;
            const auto h = std::max(hx, hy);
            const auto weight = std::sqrt(hx * hy);
            require(nx != ny && lx != ly && hx != hy && hx == hy / 2.0L &&
                        hx == lx / (nx + 1) && hy == ly / (ny + 1),
                    "rectangular geometry and spacing definitions",
                    {{"hx", hx, "== Lx/(nx+1)", lx / (nx + 1)},
                     {"hy", hy, "== Ly/(ny+1)", ly / (ny + 1)}});
            require(previous_hx == 0.0L || (previous_hx == 2.0L * hx && previous_hy == 2.0L * hy),
                    "both rectangular spacings halve exactly",
                    {{"previous_hx", previous_hx, "== 2*hx (or zero)", 2.0L * hx},
                     {"previous_hy", previous_hy, "== 2*hy (or zero)", 2.0L * hy}});

            const auto matrix = assemble_helmholtz_operator(grid, static_cast<double>(k));
            const auto sampled = manufactured_solution(grid, 1, 1);
            const auto rhs = continuous_manufactured_rhs(grid, 1, 1, static_cast<double>(k));
            const auto sx = std::sin(pi * hx / (2.0L * lx));
            const auto sy = std::sin(pi * hy / (2.0L * ly));
            const auto cx = std::cos(pi * hx / (2.0L * lx));
            const auto cy = std::cos(pi * hy / (2.0L * ly));
            const auto lambda_h = 4.0L * sx * sx / (hx * hx) + 4.0L * sy * sy / (hy * hy);
            const auto gap = lambda_h - k * k;
            const auto radius = 4.0L * cx * cx / (hx * hx) + 4.0L * cy * cy / (hy * hy) - k * k;
            // Geometry-specific bound; the square-specific spectral check is inapplicable.
            require(gap >= gap_lower_bound && gap_lower_bound > 0.0L && radius >= gap,
                    "rectangular nonresonant positive spectrum",
                    {{"gap", gap, ">= 43/9", gap_lower_bound}, {"radius", radius, ">= gap", gap}});
            const auto condition = radius / gap;
            const auto factor = (lambda - lambda_h) / gap;
            require(factor > 0.0L, "positive rectangular spatial error factor",
                    {{"factor", factor, ">", 0.0L}});

            Reference continuous(grid.unknown_count());
            Reference discrete(grid.unknown_count());
            Reference forcing(grid.unknown_count());
            Norms continuous_norm;
            Norms discrete_norm;
            Norms rhs_norm;
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 0; i < nx; ++i) {
                    const auto p = grid.linear_index(i, j);
                    const auto x = (i + 1) * lx / (nx + 1);
                    const auto y = (j + 1) * ly / (ny + 1);
                    require(p == j * nx + i && x == grid.x(i) && y == grid.y(j),
                            "rectangular indexing and independent coordinates",
                            {{"p", static_cast<long double>(p), "== j*nx+i", static_cast<long double>(j * nx + i)}});
                    const auto value = std::sin(pi * x / lx) * std::sin(pi * y / ly);
                    continuous[p] = value;
                    discrete[p] = (1.0L + factor) * value;
                    forcing[p] = (lambda - k * k) * value;
                    continuous_norm.add(continuous[p]);
                    discrete_norm.add(discrete[p]);
                    rhs_norm.add(Wide(rhs[p]));
                }
            }
            require(difference(sampled, continuous).linf <= 16.0L * epsilon &&
                        difference(rhs, forcing).linf <= 64.0L * epsilon * (lambda - k * k),
                    "independent rectangular continuous solution and RHS",
                    {{"sample Linf defect", difference(sampled, continuous).linf, "<=", 16.0L * epsilon},
                     {"RHS Linf defect", difference(rhs, forcing).linf, "<=", 64.0L * epsilon * (lambda - k * k)}});
            require(std::abs(continuous_norm.weighted(weight) - exact_continuous_l2) <= 16.0L * epsilon &&
                        std::abs(continuous_norm.linf - 1.0L) <= 16.0L * epsilon,
                    "rectangular mode norms and correct L2 weight",
                    {{"weighted mode norm", continuous_norm.weighted(weight), "== sqrt(Lx*Ly)/2 within rounding", exact_continuous_l2}});

            const auto result = gmres(matrix, rhs, options);
            checks::solver_status(result.status, profile);
            checks::iteration_limit(result.iterations, options.max_iterations);
            checks::nonzero_rhs(rhs_norm.l2, profile);
            const auto product = matrix.multiply(result.solution);
            const auto csr_residual_norm = difference(rhs, product);
            const auto csr_residual = csr_residual_norm.l2 / rhs_norm.l2;
            // Independent arithmetic over the stored CSR, separate from multiply().
            Norms wide_residual_norm;
            for (std::size_t row = 0; row < matrix.rows(); ++row) {
                Wide sum{};
                for (auto entry = matrix.row_ptr()[row]; entry < matrix.row_ptr()[row + 1]; ++entry) {
                    sum += Wide(matrix.values()[entry]) * Wide(result.solution[matrix.col_idx()[entry]]);
                }
                wide_residual_norm.add(Wide(rhs[row]) - sum);
            }
            const auto residual = wide_residual_norm.l2 / rhs_norm.l2;
            checks::true_residual(residual, options.tolerance);
            checks::true_residual(csr_residual, options.tolerance);

            const auto total = difference(result.solution, continuous);
            const auto algebraic = difference(result.solution, discrete);
            const auto spatial = difference(discrete, continuous);
            require(std::abs(spatial.linf - factor) <= 16.0L * epsilon &&
                        std::abs(spatial.weighted(weight) - exact_continuous_l2 * factor) <= 16.0L * epsilon,
                    "analytic rectangular spatial error norms",
                    {{"spatial L2(hx,hy)", spatial.weighted(weight), "== analytic norm within rounding", exact_continuous_l2 * factor},
                     {"spatial Linf", spatial.linf, "== factor within rounding", factor}});
            checks::algebraic_error(algebraic.weighted(weight), algebraic.linf,
                                    spatial.weighted(weight), spatial.linf, fraction, profile);
            const auto budget = wide_residual_norm.weighted(weight) / gap +
                64.0L * epsilon * condition * discrete_norm.weighted(weight);
            checks::error_budget(algebraic.weighted(weight), budget, spatial.weighted(weight), fraction, profile);
            checks::residual_separation(residual, spatial.l2, continuous_norm.l2, profile);
            const std::array<std::pair<std::string_view, long double>, 19> metrics{{
                {"hx", hx}, {"hy", hy}, {"h", h}, {"lambda", lambda},
                {"lambda_h", lambda_h}, {"spectral_gap", gap}, {"spectral_radius", radius},
                {"kappa2", condition}, {"spatial_error_factor", factor},
                {"true_residual", residual}, {"csr_true_residual", csr_residual},
                {"L2h", total.weighted(weight)}, {"Linf", total.linf},
                {"algebraic_L2h", algebraic.weighted(weight)}, {"algebraic_Linf", algebraic.linf},
                {"spatial_L2h", spatial.weighted(weight)}, {"spatial_Linf", spatial.linf},
                {"algebraic_L2h_bound", budget},
                {"gmres.final_residual", static_cast<long double>(result.final_residual)}
            }};
            for (const auto& [name, metric] : metrics) {
                require(std::isfinite(metric), "finite rectangular numerical metric: std::isfinite(metric)",
                        {{name, metric, "std::isfinite(metric) ==", 1.0L}});
            }

            long double order_l2 = 0.0L;
            long double order_inf = 0.0L;
            long double predicted_order = 0.0L;
            if (previous_hy != 0.0L) {
                require(factor < previous_factor, "analytic spatial error decreases",
                        {{"factor", factor, "< previous factor", previous_factor}});
                checks::decreasing_error(total.weighted(weight), previous_total_l2,
                                         total.linf, previous_total_inf, profile);
                order_l2 = std::log2(previous_total_l2 / total.weighted(weight));
                order_inf = std::log2(previous_total_inf / total.linf);
                predicted_order = std::log2(previous_factor / factor);
                checks::predicted_order(predicted_order);
                checks::observed_order(order_l2, predicted_order, order_slack, profile);
                checks::observed_order(order_inf, predicted_order, order_slack, profile);
            }
            std::cout << nx << ',' << ny << ',' << hx << ',' << hy << ',' << h << ','
                      << result.iterations << ",CONVERGED," << residual << ',' << csr_residual << ','
                      << total.weighted(weight) << ',' << total.linf << ',';
            if (previous_hy == 0.0L) {
                std::cout << "NA,NA,NA";
            } else {
                std::cout << order_l2 << ',' << order_inf << ',' << predicted_order;
            }
            std::cout << ',' << algebraic.weighted(weight) << ',' << algebraic.linf << ','
                      << spatial.weighted(weight) << ',' << spatial.linf << ',' << lambda << ','
                      << lambda_h << ',' << gap << ',' << condition << ',' << budget << '\n';
            previous_hx = hx;
            previous_hy = hy;
            previous_factor = factor;
            previous_total_l2 = total.weighted(weight);
            previous_total_inf = total.linf;
        }
    } catch (const std::exception& error) {
        std::ostringstream message;
        message << std::setprecision(17) << "nx=" << current_nx << ", ny=" << current_ny
                << ", hx=" << current_hx << ", hy=" << current_hy << " m; " << error.what();
        throw std::runtime_error(message.str());
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
