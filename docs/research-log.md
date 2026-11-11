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

### Stage 5A: FDFD + GMRES integration

**Date:** 2026-11-11 (local run date).

**Research question:** Can the existing public FDFD and complex GMRES APIs
recover algebraically manufactured fields with independently verified residual
and solution error?

**Hypothesis:** Away from discrete resonances, zero-start GMRES should recover
the manufactured vector to accuracy consistent with its true residual and the
system's spectral condition number.

**Implementation revision:** Base revision before these changes:
`0c9c6f1a1aeebe5eaff07fa53cecd1998eb7e6f7`. The measurements use that revision
plus the uncommitted Stage 5A changes to `CMakeLists.txt` and the new
`tests/helmholtz/test_integration.cpp`. The base hash alone does not identify
the new integration executable; preserve this patch with the measurements.

**Method:** TM_z scalar field with time convention `exp(-i omega t)`, real
lossless homogeneous material, homogeneous Dirichlet/PEC boundaries, and
`A_h = -Delta_h - k^2 I`. Grid dimensions count interior nodes and vectors use
`p = j*nx + i`. Each reference field comes from `manufactured_solution`.
The stored RHS is `b = A_h*u_exact`, computed with `CsrMatrix::multiply` and
cross-checked against `algebraic_manufactured_rhs`. GMRES starts from explicit
complex zeros. A fresh public CSR matvec computes `b - A_h*u_computed`, with
norms accumulated independently in long double. Forward error is measured
against `u_exact` separately. The check does not use GMRES's reported residual
as its acceptance oracle. Finite values, CONVERGED status, history size, and
global iteration limits are also checked. BREAKDOWN always fails these solves.

**Environment:** Ubuntu/WSL2, Linux x86_64, GCC 13.3.0
(`Ubuntu 13.3.0-6ubuntu2~24.04.1`), CMake 3.28.3, C++20 Debug,
`-g -Wall -Wextra -Wpedantic`, sequential `std::complex<double>` arithmetic.

**Solver parameters:** Tolerance `1e-12`, maximum iterations 1000, zero initial
guess in all cases. Domain lengths below are metres; k is in inverse metres.

| Case | nx x ny | Lx x Ly | hx, hy | k | Manufactured modes/amplitudes | Requested restart |
|---|---|---|---|---:|---|---:|
| small | 3 x 2 | 1 x 1 | 1/4, 1/3 | 3 | (1,2): 1+0.25i | 10 (capped at 6) |
| medium | 15 x 11 | 1 x 1 | 1/16, 1/12 | 8 | (1,2): 1+0.25i | 10 |
| altered | 7 x 5 | 1 x 0.75 | 1/8, 1/8 | 3 | (2,1): 1+0.25i | 10 |
| restarted | 5 x 4 | 1 x 1 | 1/6, 1/5 | 3 | (1,2): 1+0.25i; (1,1): 0.5-0.75i; (3,2): -0.25+0.5i | 2 |

**Predeclared checks and tolerance rationale:** Enumerate every discrete
Dirichlet mode independently using the analytical sine spectrum in long
double. Require `delta_h = min|lambda_pq,h - k^2| > 1 m^-2` and
`kappa_2 < 1000`. For this real symmetric matrix,
`kappa_2 = max|lambda_pq,h - k^2| / delta_h`. Require independent relative
residual <= `tau = 1e-12`. The solution-error threshold is
`2*kappa_2*(tau + 64*epsilon_double)`, motivated by the forward-error scale
`kappa_2*rho` and a conservative rounding allowance for RHS construction and
short CSR products. This is an engineering tolerance for these fixed cases,
not a rigorous floating-point certificate. The modal RHS cross-check uses
`64*epsilon_double*kappa_2`. Norm weights `hx*hy` cancel in the relative L2
error on each uniform grid.

**Measured results:** Final successful runs produced the following rounded
values; no spatial-discretization error is inferred.

| Case | nnz | delta_h (m^-2) | kappa_2 | Status | Iterations | Independent relative residual | Relative solution error |
|---|---:|---:|---:|---|---:|---:|---:|
| small | 20 | 9.372583002 | 7.748922254 | CONVERGED | 1 | 2.803805073e-16 | 1.922962686e-16 |
| medium | 773 | 13.55836306 | 111.8386264 | CONVERGED | 1 | 2.951992605e-13 | 7.156105209e-15 |
| altered | 151 | 17.89216815 | 26.60984559 | CONVERGED | 1 | 2.569974285e-15 | 4.646580477e-16 |
| restarted | 82 | 10.19532121 | 21.16703087 | CONVERGED | 64 | 8.826281448e-13 | 3.218596031e-12 |

The initial small-grid candidate used k=8 and failed the predeclared spectral
gap guard (gap < 1 m^-2). Its k was changed to 3 before the final runs; no
residual or solution-error tolerance was relaxed. The other configurations
passed unchanged. Final full CTest: 9/9 passed, with no compiler warnings
observed in the build output.

**Interpretation:** The single-mode fields are discrete eigenvectors, making
one-step convergence expected in exact arithmetic; this is a limited baseline.
The three-mode configuration required 64 iterations with restart 2, explicitly
testing continuation across cycles. Both the independent residual and reference
error satisfy the predeclared checks in these runs.

**Limitations:** These are small algebraic systems with known fields. The
reference vector need not exactly solve the stored rounded system because
forming b incurs roundoff. The independent residual still uses a double CSR
product; only its norm/difference accumulation is wider. Condition numbers are
analytical-spectrum estimates for the ideal stencil, not a numerical SVD of
the rounded matrix. This campaign does not verify Stage 5B spatial convergence,
near-resonant behavior, scalability, or a user-facing Helmholtz CLI.

**Reproduction commands:** From the repository root with this Stage 5A patch:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
ctest --test-dir build -R '^wavecore.integration' -V
```

The integration executable prints status, iterations, requested restart,
spectral gap, condition number, true residual, reference error, error threshold,
and the full residual history. Raw verbose output from the recorded run is
saved locally in `build/stage5a-results.txt` (generated artifact). Reproduce a
single case with `./build/wavecore_integration_test restarted`.

**Related publication:** No new manuscript or published result for Stage 5A.

## Known Limitations

- The current records cover diagnostic linear systems and small algebraic
  FDFD + GMRES integration tests; a user-facing application workflow is absent.
- No spatial-refinement or convergence-order measurement has been recorded.
- Near-resonance conditioning has not been systematically swept.
- No preconditioning, large-scale scalability, or performance result is
  available.
- Reported values are revision- and environment-specific observations.

## Open Research Questions

- How does the verified algebraic FDFD + GMRES baseline extend to continuous
  forcing and a user-facing reproducible workflow?
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
