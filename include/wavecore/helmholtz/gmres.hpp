#pragma once

#include "wavecore/helmholtz/csr_matrix.hpp"

#include <cstddef>
#include <vector>

namespace wavecore::helmholtz {

struct GmresOptions {
    double tolerance = 1.0e-10;
    std::size_t restart = 30;
    std::size_t max_iterations = 1000;
};

enum class GmresStatus { CONVERGED, NOT_CONVERGED, BREAKDOWN };

struct GmresResult {
    CsrMatrix::value_array solution;
    GmresStatus status = GmresStatus::NOT_CONVERGED;
    std::size_t iterations = 0;
    // ||b-Ax||_2 / ||b||_2; absolute ||b-Ax||_2 when b is exactly zero.
    double final_residual = 0.0;
    // True residual at x0, followed by one entry per Arnoldi iteration.
    std::vector<double> residual_history;
};

// Sequential, unpreconditioned restarted GMRES for square complex systems.
// An empty initial_guess selects the zero vector. restart is capped at n.
// tolerance must be finite and positive; restart and max_iterations must be positive.
// Invalid inputs throw invalid_argument; nonfinite arithmetic throws overflow_error.
// BREAKDOWN means no further Krylov direction or a singular triangular problem
// was encountered without meeting the true-residual tolerance.
//
// Numerical limits:
// - A small residual does not bound the forward solution error for ill-conditioned A.
// - Arnoldi uses two MGS passes; near-breakdown is detected when the remaining
//   direction norm is <= 16 * epsilon * ||A*v||. It may cause BREAKDOWN even for
//   nonsingular A when the requested tolerance is below attainable precision.
// - Restart can stagnate. Stagnation has no separate status: the solver continues
//   until max_iterations (NOT_CONVERGED), true convergence, or BREAKDOWN.
// - QR residual estimates never determine success; b-A*x is recomputed each step.
//   This residual is still evaluated in double precision and can suffer cancellation.
// - Finite inputs do not guarantee representable intermediates. Arithmetic
//   overflow raises overflow_error; underflow or extreme dynamic range may lose
//   information. No universal scale-invariance or forward-stability guarantee.
// - Workspace is O(n*m + m*m); true-residual checks add a sparse product each step.
[[nodiscard]] GmresResult gmres(
    const CsrMatrix& matrix,
    const CsrMatrix::value_array& rhs,
    const GmresOptions& options = {},
    const CsrMatrix::value_array& initial_guess = {});

}  // namespace wavecore::helmholtz
