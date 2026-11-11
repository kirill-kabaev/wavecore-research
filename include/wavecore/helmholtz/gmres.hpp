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
[[nodiscard]] GmresResult gmres(
    const CsrMatrix& matrix,
    const CsrMatrix::value_array& rhs,
    const GmresOptions& options = {},
    const CsrMatrix::value_array& initial_guess = {});

}  // namespace wavecore::helmholtz
