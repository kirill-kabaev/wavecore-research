# Research Log

## Purpose and Conventions

This log records reproducible WaveCore Research milestones, observations, and
open questions. It separates measured observations from interpretations and
hypotheses. A numerical value is included only when it is reported in the
GMRES manuscript and tied to the executable test revision. The log does not
claim general stability from a finite set of deterministic tests.

For a reusable record structure, see the [Numerical Experiment Template](templates/numerical-experiment.md).

The experiment revision is the solver/test revision
`26bf63a73a7201c36d6e6d6933942624ef74d9de`. It must not be replaced by a later
documentation or manuscript commit when reproducing the measurements.

## Research Timeline

The repository history records the following completed stages. No calendar
dates are inferred here.

1. **Helmholtz FDFD mathematical foundation.** The project specified the
   two-dimensional problem, homogeneous Dirichlet boundaries, manufactured
   modes, and continuous versus discrete right-hand sides.
2. **Grid, CSR, and assembly implementation.** `Grid2D`, complex CSR storage,
   and the five-point `-Delta_h - k^2 I` operator were implemented with
   deterministic unit tests.
3. **Restarted complex GMRES.** The sequential unpreconditioned solver was
   implemented with Hermitian Arnoldi products, two-pass MGS, complex Givens
   rotations, restart control, and true-residual verification.
4. **Numerical GMRES verification.** Extended deterministic diagnostics were
   added for conditioning, nonnormality, restart, breakdown, scaling, residual
   gaps, and Arnoldi orthogonality.
5. **Preliminary GMRES research paper.** The measured observations were
   documented in the manuscript *Numerical Robustness of Restarted GMRES for
   Complex-Valued Linear Systems: Preliminary Computational Experiments*.

The numerical records below refer specifically to stage 4 and the manuscript
prepared from it.

## Experiment Records

### Restarted complex GMRES diagnostic campaign

**Research question:** How do restart, conditioning, nonnormality, finite
scaling, and Arnoldi diagnostics affect a small sequential complex GMRES
implementation?

**Hypothesis:** A true-residual check should prevent an optimistic
least-squares estimate from being reported as convergence, while residual size
and forward solution error can differ substantially for ill-conditioned
systems. Restart and breakdown should be reflected in the returned status and
iteration count.

**Method:** Deterministic C++ tests exercised the public GMRES API and an
internal test observer. The production method uses the Hermitian product,
two-pass modified Gram--Schmidt, complex Givens rotations, restarted Krylov
cycles, and recomputation of `b - A*x`. Independent dense or wide-precision
calculations checked residuals, reference-vector errors, and
`||V^H V-I||_F`.

**Implementation revision:**
`26bf63a73a7201c36d6e6d6933942624ef74d9de`.

**Environment:** GCC 13.3.0, CMake 3.28.3, x86-64 Ubuntu under WSL2, C++20
Debug build with `-g -Wall -Wextra -Wpedantic`. Production arithmetic uses
double-precision `std::complex`; diagnostic accumulations use
`std::complex<long double>` where stated. Runs are sequential and deterministic;
no random seed is required.

**Solver parameters:** The principal systems use tolerance `10^-12` unless
otherwise noted. Complex Hilbert(10), near-singular diagonal, and Jordan-like
systems use restart lengths 10, 2, and 6, with maximum iterations 30, 4, and
12 respectively. The scaled 3-by-3 system uses restart 99 (capped at 3) and
maximum 6 iterations. GMRES(1) stagnation uses restart 1 and maximum 5. The
residual-gap repeat uses tolerance `10^-30`, restart 3, and maximum 6.

**Measured results:** Values below are transcribed from the manuscript's
independently recomputed tables.

| System | Iterations | Status | Relative solution error `eta_x` | Independent residual `rho` | Orthogonality defect `delta_V` |
|---|---:|---|---:|---:|---:|
| Complex Hilbert(10) | 9 | CONVERGED | `5.473e-2` | `5.610e-15` | `5.453e-16` |
| Near-singular diagonal | 1 | CONVERGED | `7.071e-1` | `1.000e-14` | `0` |
| Nonnormal Jordan-like | 6 | CONVERGED | `2.860e-15` | `1.368e-16` | `4.526e-16` |
| GMRES(1) stagnation | 5 | NOT_CONVERGED | `1` | `1` | `0` |

The Hilbert one-pass control had final nine-vector defect `1.022e-4` after
eight Arnoldi steps, whereas the production two-pass diagnostic reached maximum
defect `5.453e-16`. The scaled system converged in three iterations for each
scale:

| Scale | `eta_x` | `rho` | Estimated residual | `delta_V` |
|---:|---:|---:|---:|---:|
| `1e-200` | `4.276e-16` | `2.716e-16` | `1.534e-33` | `2.076e-16` |
| `1` | `1.711e-16` | `1.608e-16` | `4.284e-32` | `6.618e-16` |
| `1e200` | `1.763e-16` | `1.443e-16` | `2.995e-33` | `2.913e-16` |

**Interpretation:** The near-singular diagonal converged by the residual
criterion while retaining a `0.7071` forward error, demonstrating the expected
conditioning distinction. The Hilbert case likewise had a very small residual
with a `5.47%` solution error. The Jordan-like case recovered its reference
solution to the reported error scale. GMRES(1) stagnated at residual one, while
the manuscript's restart-two skew-system control converged in two iterations.
The selected finite scaling cases converged in three iterations.

**Observation versus proof:** These are observations for the listed matrices,
parameters, compiler, and revision. They are not a proof of scale invariance,
universal orthogonality preservation, or general GMRES stability.

**Limitations:** The systems are small and deterministic. Only one platform was
reported; there is no cross-platform rounding study, performance benchmark,
preconditioner, condition estimator, or large Helmholtz solve. The true
residual is recomputed in double precision and can still be affected by
cancellation, underflow, or overflow. The one-pass and two-pass orthogonality
diagnostics use different aggregation rules and are a focused comparison, not a
formal rounding-error analysis.

**Reproduction commands:** From the repository root at the stated revision:

```bash
git rev-parse HEAD
git diff 26bf63a73a7201c36d6e6d6933942624ef74d9de -- \
  CMakeLists.txt include/wavecore/helmholtz/gmres.hpp \
  src/helmholtz/gmres.cpp src/helmholtz/gmres_diagnostics.hpp \
  tests/helmholtz/test_gmres.cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
./build/wavecore_gmres_test
```

The manuscript also records GCC, CMake, platform, flags, and the diagnostic
output labels needed to interpret the run.

**Related publication:** <a href="../papers/02_gmres_numerical_observations/main.tex">`papers/02_gmres_numerical_observations/main.tex`</a>

## Known Limitations

- The current records concern diagnostic linear systems, not an integrated
  FDFD Helmholtz application.
- No spatial-refinement or convergence-order measurement has been recorded.
- Near-resonance conditioning has not been systematically swept.
- No preconditioning, large-scale scalability, or performance result is
  available.
- Reported values are revision- and environment-specific observations.

## Open Research Questions

- How does the assembled FDFD Helmholtz operator behave when solved through
  the public GMRES API as one reproducible workflow?
- What spatial convergence order is observed when discretization and iteration
  errors are separated?
- How do continuous and discrete resonance proximity affect conditioning,
  residuals, and forward error?
- Which preconditioners improve Helmholtz iteration counts without obscuring
  numerical equivalence?
- How do memory use, iteration counts, and runtime scale for larger sparse
  systems and across supported platforms?

## How to Add a New Entry

Copy this template and keep every field explicit. Use `Not recorded` when the
source does not provide a parameter or measurement.

```text
### <Short experiment title>

Research question:
Hypothesis:
Method:
Implementation revision:
Environment:
Solver parameters:
Measured results:
Interpretation:
Limitations:
Reproduction commands:
Related publication:
```

Separate measured results from interpretation, identify hypotheses as such, and
link every numerical value to a source revision and reproducible command.
