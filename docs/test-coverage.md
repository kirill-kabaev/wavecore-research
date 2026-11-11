# Numerical Core Test Coverage

This page documents the tests that are actually registered by
`CMakeLists.txt`. The matrix describes assertions in the
test sources, not an inferred percentage of code coverage. A passing test is
evidence for the exercised cases only; it is not a general stability or
correctness proof for all numerical inputs.

## Registered test matrix

| Test / CTest name | Component | Verified properties | Test methodology | Limitations |
|---|---|---|---|---|
| `wavecore.library_metadata` / `test_library.cpp` | Library metadata | Project name equals `WaveCore Research`; compiled version equals the CMake-provided `WAVECORE_EXPECTED_VERSION`. | Direct equality checks in `main`; nonzero return on mismatch. | Does not exercise numerical kernels or CLI behavior. |
| `wavecore.grid` / `test_grid.cpp` | `Grid2D` | Rectangular dimensions and physical lengths; `nx`, `ny`, `unknown_count`; `hx`, `hy`; first/last interior coordinates; x-fastest linear indexing; bijective inverse mapping; out-of-range checks; invalid dimensions/lengths; size overflow. | Deterministic `check` assertions and typed exception checks for representative valid and invalid inputs. | Does not test very large valid grids, floating-point coordinate conditioning, or performance. |
| `wavecore.csr_matrix` / `test_csr_matrix.cpp` | Complex `CsrMatrix` | Known complex square product; rectangular product; zero-by-zero matrix; empty rows; CSR size/order/bounds/value validation; vector-size validation; multiplication preserves matrix and input vector. | Explicit expected complex products, metadata checks, and `check_throws` assertions. | Does not measure sparse performance, memory use, or arithmetic accuracy across broad scales. |
| `wavecore.assembly` / `test_assembly.cpp` | Helmholtz FDFD assembly and manufactured data | Control `nnz` counts for 1x1, 2x2, 3x2, and 5x5 grids; rectangular stencil and coefficients; sorted/unique CSR columns; no row-wrap connection; transpose symmetry; continuous/discrete eigenvalues; continuous and algebraic RHS formulas; `A_h u = f_alg` on several grids; invalid parameters and overflow. | Exact small-matrix CSR comparisons, tolerance-based complex comparisons, dense symmetry reconstruction, and typed exception checks. | Uses small deterministic grids only; does not establish grid-refinement order, resonance behavior, or end-to-end solver performance. |
| `wavecore.gmres` / `test_gmres.cpp` | Restarted complex GMRES and diagnostics | Known identity, diagonal, and complex systems; Hermitian inner product behavior through solver diagnostics; two-pass MGS orthogonality; complex Givens behavior through convergence/residual checks; independently recomputed true residual; restart limits; happy breakdown; singular breakdown; zero RHS and exact initial guess; `NOT_CONVERGED`; ill-conditioned Hilbert-like system; complex nonnormal Jordan-like system; extreme finite scaling; estimated/true residual gaps; stagnation; invalid inputs and nonfinite arithmetic. | Public solver results are compared with an internal observation hook and independent wide-precision residual/error calculations. A local one-pass Arnoldi control is compared with the production two-pass basis defect. Assertions check statuses, errors, residual history, iteration limits, finite values, and diagnostics. | Small hand-constructed systems cannot characterize all conditioning, field sizes, restart choices, hardware, or long-run behavior. The diagnostic observer is a test-only detail API, not a general production instrumentation interface. |

The five unit-test names above are complemented by four integration tests:

| Test / CTest name | Component | Verified properties | Test methodology | Limitations |
|---|---|---|---|---|
| `wavecore.integration.small`, `wavecore.integration.medium`, `wavecore.integration.altered`, `wavecore.integration.restarted` / `test_integration.cpp` | End-to-end FDFD + GMRES | Nonresonant discrete spectrum; assembly and manufactured algebraic RHS consistency; zero-start solve; explicit CONVERGED status; independently recomputed CSR residual; separate relative solution error; finite solution/history; global iteration limit; actual restart in the mixed-mode case. | Compose existing public APIs. Form `b = A_h*u_exact`, cross-check the modal algebraic RHS, run GMRES, then use a fresh CSR matvec and long-double norm accumulation. Require residual <= `1e-12` and solution error <= `2*kappa_2*(1e-12 + 64*epsilon_double)`. | Small fixed grids, homogeneous Dirichlet conditions, and algebraic forcing only. Single-mode cases are eigenvectors and can converge in one step; the three-mode case exercises GMRES(2) restarts. No spatial-convergence claim. |

Three spatial-convergence experiments complement those nine registrations:

| Test / CTest name | Component | Verified properties | Test methodology | Limitations |
|---|---|---|---|---|
| `wavecore.convergence` / `test_convergence.cpp` | Continuous-RHS FDFD + GMRES grid convergence | Four exactly halved grid spacings; continuous forcing checked against an independent PDE formula; nonresonant positive spectrum; CONVERGED status; independently recomputed true residual; weighted L2 and Linf errors; decreasing errors; orders in both norms consistent with the analytic mode prediction; algebraic error and a residual/rounding-aware L2 bound below 0.1% of spatial error. | Unit square, mode (1,1), k=1, nx=ny=7,15,31,63. Reuse public assembly, manufactured data, CSR matvec and GMRES. Compare the numerical solution with both independent long-double continuous samples and the exact discrete modal solution. Require observed order in [1.99,2.03], additionally within log2(1.001/0.999) of the exact modal order. | One smooth mode, square isotropic grids, safely subresonant k, and one platform. The eigenvector RHS makes GMRES a one-step solve; this is not a solver-performance, multimode-convergence, or near-resonance study. |
| `wavecore.convergence.multimode` / `test_multimode_convergence.cpp` | Two-mode continuous-RHS FDFD + GMRES refinement | Mixture of modes (1,1) and (2,3) with amplitudes 1 and -0.3; continuous RHS independently checked; four nested square grids; positive spectral gap; CONVERGED and fresh CSR true residual; L2(h)/Linf total, spatial, and algebraic errors; discrete modal orthogonality identity for L2; actual and budgeted algebraic L2 errors below 0.1% of spatial L2 error; actual algebraic Linf below 0.1%; decreasing errors and orders consistent with separate analytic predictions; more than one GMRES iteration on at least one grid. | k=1, nx=ny=7,15,31,63; zero start, tolerance=1e-11, restart=20, max_iterations=200. Compare with independent long-double continuous samples and exact discrete modal scaling. Observed orders must lie in [1.99,2.06] and within log2(1.001/0.999) of the respective analytic orders. | Two fixed real modes on square grids, one nonresonant k, one environment. Linf maxima can move between nodes. This case does not force restart, test general multimode data, or establish solver scalability. |

The rectangular extension is registered separately:

| Test / CTest name | Component | Verified properties | Test methodology | Limitations |
|---|---|---|---|---|
| `wavecore.convergence.rectangular` / `test_rectangular_convergence.cpp` | Rectangular continuous-RHS FDFD + GMRES refinement | Lx=1, Ly=1.5, nx!=ny, hx!=hy; exact joint halving; independent indexing/coordinates, continuous RHS and discrete modal reference; full positive spectral gap; correct sqrt(hx hy) norm weight; finite metrics; CONVERGED, independent wide-CSR and fresh public-CSR residuals; separate algebraic/spatial errors; decreasing errors and second-order agreement with pre-run predictions. | Grids 7x5,15x11,31x23,63x47, mode (1,1), k=1. Reuse applicable shared acceptance functions unchanged; use a rectangle-specific analytic gap bound rather than the square-specific helper. Check spatial norm identities and observed orders against independently derived predictions with the original 0.1% contamination allowance. | One rectangle, one eigenvector RHS, fixed hx/hy=1/2, one platform. No unequal-rate refinement, near-resonance, multimode rectangular study, or solver-performance claim. |

An additional persistent acceptance regression uses the same test-only policy
functions as both convergence drivers:

| Test / CTest name | Component | Verified properties | Test methodology | Limitations |
|---|---|---|---|---|
| `wavecore.verification.acceptance` / `test_verification_checks.cpp` | Stage 5B acceptance policy | Valid metrics accepted; NOT_CONVERGED/BREAKDOWN rejected; excessive independent true residual, algebraic L2/Linf errors, error budgets, non-decreasing errors, and incompatible orders rejected; NaN/Inf rejection; inclusive endpoints and nextafter-outside boundaries; iteration limits and multiple-iteration requirement; real two-mode GMRES with max_iterations=1 cannot falsely converge. | Both existing convergence drivers call `verification_checks.hpp`. Synthetic fixtures call those same semantic checks, require the specific typed failure reason and historical condition text, and verify measurement/limit diagnostics. Unrelated exceptions propagate as failures. A separate real 7x7 continuous-RHS solve uses modes (1,1)-0.3(2,3), k=1, and a fresh CSR residual. | Synthetic metrics exercise acceptance policy, not additional PDE measurements. The limited solve is one deterministic case. This does not establish platform portability, near-resonance behavior, or universal floating-point stability. |

There are thirteen CTest registrations in total. The detailed Stage 5A parameters,
spectral gaps, measured results, tolerance rationale, and reproduction commands
are recorded in the [research log](research-log.md#stage-5a-fdfd-gmres-integration).
Stage 5B formulas, parameters, and the preserved refinement table are in the
[grid-convergence record](research-log.md#stage-5b-grid-convergence-study).
The separate [multimode record](research-log.md#stage-5b-hardening-multimode-convergence)
preserves the mixture parameters, independent predictions, measurements, and CSV path.
The [rectangular record](research-log.md#stage-5c2-rectangular-grid-convergence)
contains the predeclared protocol, analytic and measured tables, and separate CSV path.

## Component coverage

### Covered by existing tests

- **Grid2D:** construction of a non-square physical grid, spacing and
  coordinates, x-fastest `(i, j) <-> p` mapping, invalid indices, invalid
  parameters, and overflow rejection.
- **Complex CSR:** complex matrix-vector multiplication, rectangular and empty
  matrices, empty rows, immutable/read-only behavior as observed through input
  preservation, and strict structural/value validation.
- **Helmholtz FDFD assembly:** five-point coefficients, boundary-neighbour
  omission, canonical CSR ordering, expected nonzero counts, absence of false
  row-wrap links, and symmetry.
- **Manufactured solutions and RHS:** continuous and discrete eigenvalue
  formulas, complex amplitudes, distinction between continuous and algebraic
  RHS, and algebraic identity of `A_h` applied to the manufactured mode.
- **Restarted complex GMRES:** complex systems with known solutions, true
  residual confirmation, residual history shape and finiteness, restart and
  global iteration limits, happy/breakdown cases, zero-RHS handling, and
  invalid-input behavior.
- **FDFD + GMRES integration:** four deterministic configurations compose
  assembly, manufactured data, CSR multiplication, and the public solver;
  success requires both true-residual and solution-error checks.
- **Spatial convergence:** continuous forcing on four nested square grids
  verifies weighted L2 and Linf second-order errors for mode (1,1) and a
  separate mixture (1,1)-0.3(2,3), with
  independent residuals and explicit solver/discretization-error separation.
- **Rectangular convergence:** a separate four-grid study verifies continuous
  forcing, physical-length scaling, unequal coordinate spacings, weighted norms,
  a geometry-specific spectral bound, and spatial order for mode (1,1).
- **Persistent acceptance regressions:** single-mode and multimode policy
  rejection paths are registered in CTest, including finite-value checks and
  boundary fixtures. These replace reliance on a one-off GDB audit for those
  paths; the historical Stage 5B CSV files and acceptance thresholds are unchanged.
  Diagnostic regressions require the complete `condition` field before
  `; measured[`, reject prefix-only matches and missing delimiters, and retain
  the typed failure reason as the primary contract. Helper preconditions and
  first-failure ordering are documented in the test-only header.

### Partially covered

- **GMRES Hermitian inner product:** complex diagnostic systems exercise the
  conjugated first argument and verify the resulting solver behavior, but no
  standalone public inner-product test exists.
- **Two-pass MGS:** the diagnostic observer computes `||V^H V-I||` for the
  production basis and compares it with a local one-pass control. This is a
  focused diagnostic, not a broad floating-point study.
- **Complex Givens rotations:** complex convergence and residual assertions
  exercise the implementation indirectly; the private rotation helper is not
  exposed or tested through a separate public unit-test API.
- **True versus estimated residual:** the tests deliberately find an optimistic
  internal estimate and assert that it cannot produce false `CONVERGED`, but
  only a small set of systems is used.
- **Ill-conditioned, nonnormal, and scaled systems:** Hilbert-like,
  Jordan-like, and `1e-200`/`1e200` scale cases are covered, without a broad
  parameter sweep or cross-platform comparison.
- **FDFD/GMRES relationship:** algebraic manufactured systems now have dedicated
  integration tests, including a mixed-mode restarted solve. Continuous forcing
  and spatial refinement cover square single-mode/two-mode fields and one
  rectangular single-mode case. Broader mixtures/geometries, unequal-rate
  directional refinement, and an application-level workflow remain gaps.

### Not yet covered

- Spatial convergence for broader multimode data, other rectangular aspect
  ratios, unequal-rate directional refinement, and indefinite Helmholtz operators.
- Near-resonant Helmholtz parameters and systematic eigenvalue proximity tests.
- Preconditioned solver behavior (no preconditioner exists yet).
- Large sparse matrices, memory pressure, and scalability.
- Cross-platform/compiler reproducibility.
- Performance, throughput, and iteration-time regressions.
- Long-running or randomized property-based numerical campaigns.

## Future Test Coverage

The following gaps require new tests rather than stronger interpretation of the
current assertions:

| Gap | Future test type |
|---|---|
| Extended FDFD+GMRES integration | Extend the algebraic, single-mode, and two-mode continuous-RHS baselines to an application-level input/output workflow, preserving independent residual and solution-error checks. |
| Broader grid refinement and convergence order | Extend the square single-mode/two-mode and rectangular single-mode studies to broader mixtures, other aspect ratios, unequal-rate directional refinement, and indefinite nonresonant operators with conditioning-aware solver tolerances. |
| Near-resonance behavior | Sweep `k^2` around continuous and discrete eigenvalues, recording conditioning indicators, solver status, and residual/error separately. |
| Preconditioning | Add controlled solver tests comparing unpreconditioned and approved preconditioned variants on identical matrices and stopping criteria. |
| Large sparse systems | Generate reproducible larger CSR systems and test allocation limits, memory use, iteration counts, and failure handling. |
| Cross-platform verification | Repeat deterministic builds/tests with supported GCC and Clang configurations and compare toleranced numerical outputs. |
| Performance and scalability | Add benchmark-style measurements with fixed compiler, hardware, problem size, and thread metadata; keep them separate from correctness tests. |

No numerical percentages or benchmark results are asserted here because the
repository currently provides no coverage report or performance measurement
artifact for those quantities.
