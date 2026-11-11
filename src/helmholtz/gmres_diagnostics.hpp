#pragma once

// Internal test instrumentation. Not part of the installed/public solver API.
#include "wavecore/helmholtz/gmres.hpp"

#include <functional>

namespace wavecore::helmholtz::detail {

// Called after each nonsingular least-squares solve, before convergence checks.
// basis contains the actual Arnoldi vectors used by this iteration; next contains
// the unnormalized new direction and is excluded if breakdown was detected.
using ArnoldiObserver = std::function<void(
    const std::vector<CsrMatrix::value_array>& basis,
    const CsrMatrix::value_array& next,
    double next_norm,
    bool breakdown,
    double estimated_residual,
    double true_residual)>;

GmresResult gmres_observed(const CsrMatrix& matrix,
                           const CsrMatrix::value_array& rhs,
                           const GmresOptions& options,
                           const CsrMatrix::value_array& initial_guess,
                           const ArnoldiObserver& observer);

}  // namespace wavecore::helmholtz::detail
