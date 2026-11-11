# Research Log

## Purpose and Conventions

This log records reproducible WaveCore Research milestones, observations, and
open questions. It separates measured observations from interpretations and
hypotheses. A numerical value is included only when it is tied to an
executable test revision and recorded output or the GMRES manuscript. The log does not
claim general stability from a finite set of deterministic tests.

For a reusable record structure, see the [Numerical Experiment Template](templates/numerical-experiment.md).

The original GMRES diagnostic campaign uses the solver/test revision
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
6. **Stage 5A integration.** Four deterministic algebraic FDFD + GMRES cases
   independently verify residuals and solution errors.
7. **Stage 5B grid convergence.** Continuous forcing on four nested grids
   confirms second-order L2/Linf accuracy for one smooth nonresonant mode,
   with separate algebraic and discretization errors.
8. **Stage 5B multimode hardening.** An independent two-mode experiment on the
   same grids checks both norm orders against modal predictions and exercises
   more than one GMRES iteration. The original single-mode experiment is preserved.

The numerical records below cover the GMRES campaign and Stages 5A–5B;
each record identifies its own implementation revision and scope.

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

### Stage 5B: Grid convergence study

**Date:** 2026-11-11 (local run date).

**Research question:** Does the five-point FDFD scheme exhibit second-order
spatial accuracy when the continuous RHS is used and algebraic error is controlled?

**Hypothesis:** For a fixed smooth nonresonant Dirichlet mode, both weighted
L2 and Linf errors decay as h squared. This hypothesis concerns spatial accuracy,
not general GMRES robustness or efficiency.

**Implementation revision:** Base
`69c27752fc022a04d200fa4f176dc9011db65b59` plus the uncommitted Stage 5B
`CMakeLists.txt` registration and `tests/helmholtz/test_convergence.cpp`.
The numerical-core implementation is unchanged from that base.

**Environment:** Ubuntu 24.04 under WSL2, Linux
`6.18.58.1-microsoft-standard-WSL2`, x86_64; GCC 13.3.0; CMake 3.28.3;
C++20 Debug, `-g -std=c++20 -Wall -Wextra -Wpedantic`. Assembly/GMRES use
`std::complex<double>`; independent reference formulas and norm accumulation use
`std::complex<long double>` / `long double`. Sequential, no preconditioner or
external numerical libraries. No performance benchmark is claimed.

**Mathematical model and parameters:** Homogeneous isotropic lossless TM_z,
time convention exp(-i omega t), field u=E_z, unit square in metres, zero
Dirichlet/PEC boundaries. This study intentionally uses a different mode and
wavenumber from the initial (1,2), k=8 specification example:

$$
-\Delta u-k^2u=f,\qquad k=1\ \mathrm{m}^{-1},\qquad
u^*(x,y)=\sin(\pi x/L_x)\sin(\pi y/L_y),\quad L_x=L_y=1\ \mathrm m.
$$

Amplitude is 1 V/m; f has units V/m cubed. The existing
`continuous_manufactured_rhs(grid,1,1,k)` supplies
f=(2 pi squared - k squared) u*, with the eigenvalue in inverse square metres.
It is independently checked against the PDE formula in long double; **b is
not formed as A_h*u_exact**. Grid dimensions count interior nodes:
nx=ny=7,15,31,63, h=1/(n+1)=1/8,1/16,1/32,1/64 m. Coordinates are
x=(i+1)hx, y=(j+1)hy, p=j*nx+i. Boundaries are omitted from the unknown vector.

**Nonresonance and conditioning:** The continuous minimum eigenvalue is
2 pi squared, so the continuous gap is 2 pi squared - 1 > 0. On these square grids,

$$
\lambda_{11,h}=\frac{8}{h^2}\sin^2\!\left(\frac{\pi h}{2}\right),\qquad
\delta_h=\lambda_{11,h}-1,\qquad
\kappa_2(A_h)=\frac{8h^{-2}\cos^2(\pi h/2)-1}{\delta_h}.
$$

Here and below lengths are numerical values in metres. The minimum/maximum
discrete eigenvalues occur at (1,1)/(n,n). Since sin(z)>=2z/pi on [0,pi/2],
lambda_11,h>=8 and delta_h>=7 inverse square metres. Thus **all** eigenvalues
are positive and separated from zero, not just the forced mode. The test checks
this bound using independent formulas rather than the production eigenvalue helper.
Measured gaps range from 18.4868396771 to 18.7352455345; condition numbers grow
from 26.5872 to 1747.8962. This is a subresonant positive-definite example.

**Method and error definitions:** Run the public CSR assembly and GMRES;
recompute b-A_h*u_computed with a fresh public CSR matvec, accumulating norms
independently. The relative true residual is ||b-A_h*u_computed||_2/||b||_2;
the weighted-norm ratio is identical on each uniform grid. For total nodal
error e=u_computed-u*|h,

$$
E_{2,h}=\left(h_xh_y\sum_p|e_p|^2\right)^{1/2},\qquad
E_\infty=\max_p|e_p|,\qquad
p_\ell=\log_2(E_{\ell-1}/E_\ell).
$$

These are **absolute** errors: E_2,h has units V, E_inf has units V/m.
Orders are computed separately for each norm. In exact arithmetic the discrete
solution and its spatial-error factor are

$$
u_h=\frac{2\pi^2-1}{\lambda_{11,h}-1}u^*|_h=(1+d_h)u^*|_h,\qquad
d_h=\frac{2\pi^2-\lambda_{11,h}}{\lambda_{11,h}-1}>0.
$$

Independent long-double samples of u* and u_h separate the measured algebraic
error u_computed-u_h from spatial error u_h-u*. On these odd-n grids the centre
is sampled, ||u*||_inf=1 and ||u*||_2,h=1/2; hence the exact spatial errors are
d_h and d_h/2. The test checks these identities. Stored-matrix and RHS
rounding are included in the numerical/algebraic error budget; the reference
u_h is the exact-arithmetic stencil solution rather than an exact solution of
every rounded input. For these binary h and k the stencil coefficients themselves
are exactly representable.

**Solver parameters and acceptance criteria:** Zero initial guess;
unpreconditioned restarted GMRES, restart=20, max_iterations=200,
tolerance=1e-11. Success requires CONVERGED, finite values, the iteration limit,
and independently computed true residual <=1e-11. This residual threshold is
accepted only with the following forward-error checks:

- Actual algebraic errors in **both** norms must be <=0.001 times the respective
  exact spatial errors. This limits contamination of the measured orders.
- A conditioning-aware L2 budget must also be <=0.001 times spatial L2 error:
  B_2=||r||_2,h/delta_h + 64 epsilon_double kappa_2 ||u_h||_2,h.
  Actual algebraic L2 error must be <=B_2. The first term follows from
  ||A_h inverse||_2=1/delta_h; the second is a conservative allowance for
  constructing the references/forcing and rounding sparse products, not a
  universal floating-point error theorem.
- True relative residual must be at least 100 times below the relative spatial
  error, and both measured total errors must strictly decrease on every refinement.

**Order rationale:** Taylor expansion gives

$$
\lambda_{11,h}=2\pi^2-\frac{\pi^4}{6}h^2+\frac{\pi^6}{180}h^4+O(h^6),
\qquad
d_h=\frac{\pi^4}{6(2\pi^2-1)}h^2[1+O(h^2)].
$$

The exact modal orders $\log_2(d_h/d_{h/2})$, evaluated from the trigonometric
formula, are 2.0091162938184, 2.0022732186790, 2.0005679399975. The test checks
that these lie in [2,2.02]. At <=0.1% algebraic contamination per level, the
order perturbation is bounded by log2(1.001/0.999), approximately 0.002886.
Observed orders must satisfy that tighter distance to their analytic prediction
as well as [1.99,2.03]. These bounds follow from this mode and grid sequence,
not from an arbitrary generic second-order threshold. They were set before
the first run; no tolerance was loosened to obtain a pass.

**Measured results:** All four solves returned CONVERGED. Errors below are
total numerical-minus-continuous errors, not algebraic-only errors.

| Interior grid | h (m) | Iterations | True relative residual | L2(h) error (V) | Linf error (V/m) | Order L2 / Linf |
|---|---|---|---|---|---|---|
| 7 x 7 | 0.125 | 1 | 1.7385039704e-15 | 6.8256427133e-3 | 1.3651285427e-2 | N/A |
| 15 x 15 | 0.0625 | 1 | 7.3024031066e-15 | 1.6956619789e-3 | 3.3913239578e-3 | 2.0091162938 / 2.0091162938 |
| 31 x 31 | 0.03125 | 1 | 3.2314534745e-14 | 4.2324806759e-4 | 8.4649613518e-4 | 2.0022732187 / 2.0022732187 |
| 63 x 63 | 0.015625 | 1 | 1.3350198527e-13 | 1.0577037050e-4 | 2.1154074100e-4 | 2.0005679400 / 2.0005679400 |

| Interior grid | Algebraic L2(h) | Algebraic Linf | L2 error budget B_2 |
|---|---|---|---|
| 7 x 7 | 9.5100590592e-17 | 2.1703016988e-16 | 1.9237339170e-13 |
| 15 x 15 | 3.0962587702e-16 | 6.8651681562e-16 | 7.7759562209e-13 |
| 31 x 31 | 2.1030318917e-16 | 5.5662939535e-16 | 3.1197571542e-12 |
| 63 x 63 | 4.7643231173e-16 | 1.1758714662e-15 | 1.2488941937e-11 |

Full-precision rows, including both orders, analytic predictions, separate
errors, spectral gaps, and condition numbers, are preserved in the repository:
`experiments/helmholtz_2d/results/stage5b_grid_convergence.csv`.
The CSV is raw executable output after its first metadata line; table entries
above are rounded for display, not additional measurements.

**Interpretation:** The measured orders approach two and agree with exact modal
predictions in both norms. Spatial errors dominate both measured algebraic
errors and their conservative L2 budgets, confirming the expected spatial order
for this experiment. One GMRES step is expected: continuous forcing is a scalar
multiple of a discrete eigenvector. It is not evidence for constant iteration
counts on general Helmholtz systems.

**Limitations:** One mode, square isotropic grids, k=1, homogeneous medium, PEC
walls, double precision, and one GCC/WSL2 environment. No indefinite/near-resonant
or multimode spatial study, preconditioning, performance/scalability conclusion,
or mesh-independent GMRES claim. Very fine grids may raise cancellation and
conditioning effects. Stage 5A retains a separate mixed-mode restart test.

**Reproduction commands:** From the repository root with the Stage 5B patch:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
ctest --test-dir build -R '^wavecore.convergence$' -V
./build/wavecore_convergence_test
```

The direct executable prints a versioned metadata line and the full CSV table;
nonzero exit status signals any acceptance failure. To keep a local raw capture:

```bash
set -o pipefail
./build/wavecore_convergence_test | tee build/stage5b-grid-convergence.txt
```

**Related publication:** No article has been modified; this is a reproducible
research-log experiment, not a peer-reviewed result.

### Stage 5B hardening: Multimode convergence

**Date:** 2026-11-11 (local run date).

**Research question and hypothesis:** Does second-order spatial convergence
persist for a smooth mixture of two distinct eigenmodes, when GMRES must perform
multiple iterations? The hypothesis is O(h squared) error in both norms after
separating numerical/algebraic error from exact-arithmetic spatial error.

**Implementation revision:** Base HEAD
`69c27752fc022a04d200fa4f176dc9011db65b59`, the existing uncommitted Stage 5B
patch, and this hardening patch (`test_multimode_convergence.cpp` and its CMake
registration). Production FDFD/GMRES is unchanged. The single-mode test and CSV
remain byte-for-byte unchanged from the preceding Stage 5B experiment.

**Environment:** Same GCC 13.3.0 / CMake 3.28.3, Ubuntu/WSL2 x86_64 C++20
Debug environment as the single-mode record. Production arithmetic is complex
double; independent reference formulas and norm accumulation use long double.
Build flags remain `-g -std=c++20 -Wall -Wextra -Wpedantic`. Sequential,
unpreconditioned, no external numerical library; no timing claim.

**Model and parameters:** Unit square in metres, TM_z field E_z with
exp(-i omega t), homogeneous isotropic lossless medium, zero Dirichlet/PEC,
A_h=-Delta_h-k squared I, k=1 inverse metre. Define
phi_mn=sin(m pi x/Lx) sin(n pi y/Ly), Lx=Ly=1 m. Amplitudes are 1 and -0.3 V/m:

$$
u^*=\phi_{11}-0.3\phi_{23},\qquad
f=(2\pi^2-1)\phi_{11}-0.3(13\pi^2-1)\phi_{23}.
$$

Eigenvalue coefficients above are numerical values in inverse square metres;
f has units V/m cubed. Existing `manufactured_solution` and
`continuous_manufactured_rhs` generate each mode, then their vectors are summed.
No forcing is constructed as A_h*u_exact. Independent coordinate-based
long-double evaluation checks both public helper results against the PDE.
Interior counts nx=ny=7,15,31,63, h=1/8,1/16,1/32,1/64 m; x=(i+1)h,
y=(j+1)h, p=j*nx+i. No random inputs.

**Independent discrete reference:** With numerical lengths in metres,

$$
\lambda_{mn}=(m^2+n^2)\pi^2,\qquad
\lambda_{mn,h}=\frac4{h^2}\left[\sin^2\frac{m\pi h}{2}
                                  +\sin^2\frac{n\pi h}{2}\right],\qquad
d_{mn,h}=\frac{\lambda_{mn}-\lambda_{mn,h}}{\lambda_{mn,h}-1},
$$

$$
u_h=(1+d_{11,h})\phi_{11}-0.3(1+d_{23,h})\phi_{23}.
$$

References use independent trigonometric formulas, not the production discrete
eigenvalue helper. The full operator is positive definite: the same sin(z)>=2z/pi
argument as in the single-mode record gives delta_h=lambda_11,h-1>=7 inverse
square metres. The two forced eigenvalues are distinct. kappa_2 is computed from
the extreme (1,1)/(n,n) modes and ranges from 26.5872 to 1747.8962.

**Norms and stopping criteria:** Total error is u_computed-u*, algebraic error
is u_computed-u_h, and spatial error is u_h-u*. Use absolute
L2(h)=sqrt(h squared sum |e| squared) and Linf=max |e|, with units V and V/m.
Discrete sine orthogonality independently gives spatial
L2(h)=sqrt(d_11,h squared + 0.09 d_23,h squared)/2; the test checks this identity.
Spatial Linf is evaluated as the maximum of the independent modal error at
the grid nodes, not a sum of individual maxima.

Zero initial guess, tolerance=1e-11, restart=20, max_iterations=200, matching
the single-mode settings. Require CONVERGED, finite diagnostics, and the global
iteration limit. Recompute true relative residual ||b-A_h*u_computed||_2/||b||_2
using a fresh public CSR product and wide norms. Require it <=1e-11 and at least
100 times below the relative spatial L2 error. At least one grid must require
more than one GMRES iteration.

The residual threshold is only accepted together with **both** actual algebraic
norm errors <=0.1% of their spatial counterparts. Also check algebraic
L2(h)<=B_2<=0.001 times spatial L2(h), where
B_2=||r||_2,h/delta_h + 64 epsilon_double kappa_2 ||u_h||_2,h.
The inverse-gap part follows from the spectrum; the rounding part is the same
conservative budget as the single-mode test, not a universal stability theorem.

**Order criteria and rationale:** For each mode,
lambda_mn-lambda_mn,h=pi to the fourth (m to the fourth + n to the fourth)
h squared/12 + O(h to the fourth), so the mixture has a nonzero second-order
leading error. Before solver measurements, independent modal calculations gave
orders approximately (2.051388,2.047707), (2.012614,1.999341),
(2.003139,2.002857) in (L2,Linf). The higher-frequency mode gives larger
coarse-grid corrections, and nodal maxima can shift, so Linf order need not
approach two monotonically or from above.

For each refinement require decreasing total errors and observed orders in
[1.99,2.06]. More tightly, each observed order must differ from its own analytic
spatial-error prediction by at most log2(1.001/0.999), approximately 0.002886.
This bound follows from the <=0.1% algebraic-error condition at both levels.
The broad interval encloses the independently calculated finite-grid orders
with this allowance; the single-mode test and its narrower interval are unchanged.
No acceptance criterion was loosened after execution.

**Measured results:** All four systems returned CONVERGED:

| Interior grid | h (m) | Iterations | True relative residual | L2(h) (V) | Linf (V/m) | Order L2 / Linf |
|---|---|---|---|---|---|---|
| 7 x 7 | 0.125 | 2 | 4.3080550329e-15 | 1.6772829107e-2 | 4.0295266132e-2 | N/A |
| 15 x 15 | 0.0625 | 2 | 9.5409393249e-14 | 4.0464757927e-3 | 9.7461426998e-3 | 2.0513881818 / 2.0477071104 |
| 31 x 31 | 0.03125 | 2 | 1.7138177559e-12 | 1.0028128572e-3 | 2.4376485621e-3 | 2.0126135665 / 1.9993411998 |
| 63 x 63 | 0.015625 | 3 | 9.1598599610e-12 | 2.5015834182e-4 | 6.0820664135e-4 | 2.0031389318 / 2.0028566712 |

| Interior grid | Algebraic L2(h) | Algebraic Linf | L2 budget B_2 |
|---|---|---|---|
| 7 x 7 | 3.3495750968e-16 | 1.0730891002e-15 | 2.0637891423e-13 |
| 15 x 15 | 1.8466325247e-15 | 5.5115959539e-15 | 9.1809021736e-13 |
| 31 x 31 | 8.9020525553e-15 | 3.4063822997e-14 | 5.1886078477e-12 |
| 63 x 63 | 1.4006429344e-14 | 6.1869712370e-14 | 2.3369923369e-11 |

Raw full-precision CSV rows, including separate analytic order predictions,
errors, spectral gaps, condition numbers, and budgets, are preserved separately:
`experiments/helmholtz_2d/results/stage5b_multimode_convergence.csv`.
This is the direct executable output after its first metadata line. The previous
`stage5b_grid_convergence.csv` remains unchanged.

**Interpretation:** Both norms show the predicted second-order spatial accuracy
for this mixture. Two distinct excited eigenvalues give a two-dimensional
Krylov space in exact arithmetic. Observed iteration counts are 2,2,2,3; finite
precision may require more than two iterations. The third step on the finest
grid is consistent with roundoff-excited components outside the ideal modal
subspace; their cause was not separately instrumented. Restart=20 is configured
but no restart is exercised here. Spatial errors dominate both measured
algebraic errors and their conservative L2 budgets.

**Limitations:** Two fixed real modes, four square isotropic grids, k=1 safely
below resonance, one compiler/platform. No general multimode, anisotropic,
indefinite, near-resonance, scalability, or universal-stability claim. No new
article or publication claim.

**Reproduction commands:** From the root with the complete Stage 5B/hardening patch:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 2
ctest --test-dir build -R '^wavecore.convergence.multimode$' -V
./build/wavecore_multimode_convergence_test
ctest --test-dir build --output-on-failure
```

### Stage 5B finalization: Diagnostic context and documentation

The finalization patch changes only failure reporting in the two convergence
drivers: errors include n, h, the violated condition, measurements, and allowed
limits. Mathematical formulas, solver parameters, and acceptance predicates
are preserved. Both original CSV files and their recorded numerical values
remain unchanged. The historical byte-for-byte statement in the multimode
record describes that earlier hardening step; finalization intentionally updates
single-mode diagnostic source code as well.

The root changelog records Stage 5A and the unreleased Stage 5B milestone.
[Technical decisions](technical-decisions.md) document the rationale evidenced
by the implemented experiments. README and navigation now reflect the completed
test workflows. No production numerical kernel or manuscript is changed.

## Known Limitations

- The current records cover diagnostic linear systems, small algebraic
  FDFD + GMRES integration tests, and single-mode/two-mode continuous-RHS refinement
  study; a user-facing application workflow is absent.
- Spatial order has been measured for one mode and one two-mode mixture on
  square nonresonant grids; broader spatial studies remain open.
- Near-resonance conditioning has not been systematically swept.
- No preconditioning, large-scale scalability, or performance result is
  available.
- Reported values are revision- and environment-specific observations.

## Open Research Questions

- How do the verified algebraic, single-mode, and two-mode continuous-RHS baselines extend
  to a user-facing reproducible workflow?
- Does the observed spatial order persist for broader modal mixtures, rectangular grids,
  anisotropic refinement, and indefinite nonresonant operators?
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
