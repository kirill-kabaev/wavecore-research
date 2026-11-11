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

The five CTest names above are the complete registrations in the current
CMake configuration.

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
- **FDFD/GMRES relationship:** each subsystem is tested independently and the
  public APIs can be composed manually, but there is no dedicated integration
  test for a complete Helmholtz solve.

### Not yet covered

- Spatial grid refinement and measured convergence order.
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
| FDFD+GMRES integration | Build a manufactured Helmholtz system, solve it through the public assembly and GMRES APIs, and independently verify the solution and true residual. |
| Grid refinement and convergence order | Run a fixed manufactured mode over a documented sequence of `nx`, `ny` grids and compute norm errors and observed orders. |
| Near-resonance behavior | Sweep `k^2` around continuous and discrete eigenvalues, recording conditioning indicators, solver status, and residual/error separately. |
| Preconditioning | Add controlled solver tests comparing unpreconditioned and approved preconditioned variants on identical matrices and stopping criteria. |
| Large sparse systems | Generate reproducible larger CSR systems and test allocation limits, memory use, iteration counts, and failure handling. |
| Cross-platform verification | Repeat deterministic builds/tests with supported GCC and Clang configurations and compare toleranced numerical outputs. |
| Performance and scalability | Add benchmark-style measurements with fixed compiler, hardware, problem size, and thread metadata; keep them separate from correctness tests. |

No numerical percentages or benchmark results are asserted here because the
repository currently provides no coverage report or performance measurement
artifact for those quantities.
